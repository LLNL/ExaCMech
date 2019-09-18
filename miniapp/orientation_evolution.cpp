#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#if defined(RAJA_ENABLE_CUDA)
#include "cuda.h"
#endif
#include "miniapp_util.h"
#include "retrieve_kernels.h"
#include "setup_kernels.h"
#include "material_kernels.h"

#include <math.h>
#include <fstream>
#include <random>
#include <sstream>
#include <string>


using namespace ecmech;

int main(int argc, char *argv[]){
   //TODO:
   //Read in options (orientation file, material model,
   //                 material model params file, device type (CPU, OpenMP, CUDA, HIP))
   //Read in data file for the orientations
   //Read in data file for material model params
   //Create material model
   //Setup all of the various variables being used on the device:
   //      state variables, stress field, all of the temporary variables,
   //      material tangent stiffness matrix, and velocity gradient.
   //
   //The below initialization should work just fine for setting up the material model
   //Place various kernels in the areas they need to be.
   //Put a giant loop around the inner kernels to allow things to evolve in time.
   //Compare first the OpenMP versus the serial implementation results
   //Compare GPU versus the serial results (We need to figure out what the bounds on our
   //differences are between our results. We should see some due to CPU and GPU being different
   //devices with potentially subtle different ways of dealing with floating point arrithmetic)

   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] <<
         " option file path which contains: quat file path, material model, "
                << "material param file path, and device type each on their own line."
                << std::endl;
      return 1;
   }

   const double dt = 0.00025;
   const int nsteps = 60;

   //All of the varibles that we'll be using in our simulations
   double* state_vars = nullptr;
   double* vgrad = nullptr;

   int nqpts = 0;
   int num_props = 0;
   int num_hardness = 0;
   int num_gdot = 0;
   int iHistLbGdot = 0;
   //We have 23 state variables plus the 4 from quaternions for
   //a total of 27 for FCC materials using either the
   //voce or mts model.
   //They are in order:
   //dp_eff(1), eq_pl_strain(2), n_evals(3), dev. elastic strain(4-8),
   //quats(9-12), h(13), gdot(14-25), rel_vol(26), int_eng(27)
   int num_state_vars = 27;

   ecmech::matModelBase* mat_model_base;
   //Could probably do this in a smarter way where we don't create two class objects for
   //our different use cases...

   ecmech::matModelEvptn_FCC_A* mat_modela = memoryManager::allocate<ecmech::matModelEvptn_FCC_A>(1);
   ecmech::matModelEvptn_FCC_B* mat_modelb = memoryManager::allocate<ecmech::matModelEvptn_FCC_B>(1);

   Accelerator device;

   //Data structures needed for each time step
   //We really don't need to allocate these constantly, so we should just do it
   //once and be done with it.
   double* stress_array = nullptr;
   double* stress_svec_p_array = nullptr;
   double* d_svec_p_array = nullptr;
   double* w_vec_array = nullptr;
   double* ddsdde_array = nullptr;
   double* vol_ratio_array = nullptr;
   double* eng_int_array = nullptr;

   #if defined(RAJA_ENABLE_CUDA)
   cudaDeviceSetLimit(cudaLimitStackSize, 16000);
   #endif

   std::string mat_model_str;

   //The below scope of work sets up everything that we're going to be doing initially.
   //We're currently doing this to ensure that memory used during the set-up is freed
   //early on before we start doing all of our computations. We could probably keep everything
   //in scope without running into memory issues.
   {
      //All the input arguments
      std::string option_file(argv[1]);

      std::string ori_file;
      std::string mat_prop_file;
      std::string device_type;

      {
         std::ifstream ofile(option_file);
         ofile.clear();
         std::string line;

         std::getline(ofile, ori_file);
         std::getline(ofile, mat_model_str);
         std::getline(ofile, mat_prop_file);
         std::getline(ofile, device_type);

         if (ofile.fail()) {
            std::cerr << "Option file could not be correctly parsed.\n"
                      << "Option file contains: quat file path, material model, "
                      << "material param file path, and device type each on their own line."
                      << std::endl;
            return 1;
         }
      }

      //Quaternion and the number of quaternions total.
      std::vector<double> quats;
      //This next chunk reads in all of the quaternions and pushes them to a vector.
      //It will exit if 4 values are not read on a line.
      {
         std::ifstream qfile(ori_file);
         std::string line;
         while (std::getline(qfile, line)) {
            std::istringstream iss(line);
            double q1, q2, q3, q4;

            nqpts += 1;

            if (!(iss >> q1 >> q2 >> q3 >> q4)) {
               std::cerr << "Quat file has malformed line on line: " << nqpts << std::endl;
               return 1;
            }    // error

            quats.push_back(q1); quats.push_back(q2); quats.push_back(q3); quats.push_back(q4);
         }
      }


      //We now detect which device is desired to run the different cases.
      //Compiler flags passed in will tell us which options are available
      //based on what RAJA was built with. If we do not have support for
      //the chosen value then we should error out and let the user know
      //which values are available.

      if (device_type.compare("CPU") == 0) {
         device = Accelerator::CPU;
      }
      #if defined(RAJA_ENABLE_OPENMP)
      else if (device_type.compare("OpenMP") == 0) {
         device = Accelerator::OPENMP;
      }
      #endif
      #if defined(RAJA_ENABLE_CUDA)
      else if (device_type.compare("CUDA") == 0) {
         device = Accelerator::CUDA;
      }
      #endif
      else {
         std::cerr << "Accelerator is not supported or RAJA was not built with" << std::endl;
         return 1;
      }

      std::cout << "About to initialize class" << std::endl;
      //Initialize our base class using the appropriate model
      if (mat_model_str.compare("voce") == 0) {
         num_props = 17;
         num_hardness = mat_modela->nH;
         num_gdot = mat_modela->nslip;
         iHistLbGdot = mat_modela->iHistLbGdot;
//For the time being I haven't included the other model yet up on the device just because I want
//to see if things will even work when using the simpler model.
#if defined(RAJA_ENABLE_CUDA)
         if (device == Accelerator::CUDA) {
            RAJA::forall<RAJA::cuda_exec<1> >(RAJA::RangeSegment(0, 1), [ = ] RAJA_HOST_DEVICE(int) {
               new(mat_modela) ecmech::matModelEvptn_FCC_A();
            });
         }
#endif
         mat_model_base = dynamic_cast<matModelBase*>(mat_modela);
      }
      else if (mat_model_str.compare("mts") == 0) {
         num_props = 24;
         num_hardness = mat_modelb->nH;
         num_gdot = mat_modelb->nslip;
         iHistLbGdot = mat_modelb->iHistLbGdot;
         mat_model_base = dynamic_cast<matModelBase*>(mat_modelb);
      }
      else {
         std::cerr << "material model must be either voce or mts " << std::endl;
         return 1;
      }

      std::cout << "Initialized class" << std::endl;

      //Read and store our material property data
      std::vector<double> mat_props;
      {
         std::ifstream mfile(mat_prop_file);
         std::string line;
         int nlines = 0;
         while (std::getline(mfile, line)) {
            std::istringstream iss(line);
            double prop;

            nlines += 1;
            // std::cout << line << std::endl;
            if (!(iss >> prop)) {
               std::cerr << "Material prop file has a malformed line on line: " << nlines << std::endl;
               return 1;
            }    // error

            mat_props.push_back(prop);
         }
         if (nlines != num_props) {
            std::cerr << "Material prop file should have " << num_props
                      << " properties (each on their own line). A total of " << nlines
                      << " properties were provided instead." << std::endl;
            return 1;
         }
      }


      //Below initializes our material model class
      //Opts and strs are just empty vectors of int and strings
      std::vector<double> params;
      std::vector<int> opts;
      std::vector<std::string> strs;

      for (int i = 0; i < mat_props.size(); i++) {
         params.push_back(mat_props.at(i));
      }

      //We really shouldn't see this change over time at least for our applications.
      mat_model_base->initFromParams(opts, params, strs);
      //
      mat_model_base->complete();

      std::cout << "Class has been completely initialized" << std::endl;

      //We're now initializing our state variables and vgrad to be used in other parts
      //of the simulations.
      state_vars = memoryManager::allocate<double>(num_state_vars * nqpts);
      vgrad = memoryManager::allocate<double>(nqpts * ecmech::ndim * ecmech::ndim);

      double* quats_array = quats.data();

      init_data(device, quats_array, mat_model_base, nqpts, num_hardness,
                num_gdot, iHistLbGdot, num_state_vars, state_vars);
      std::cout << "Data is now initialized" << std::endl;
      setup_vgrad(vgrad, nqpts);
   }

   //The stress array is the only one of the below variables that needs to be
   //initialized to 0.
   stress_array = memoryManager::allocate<double>(nqpts * ecmech::nsvec);
   for (int i = 0; i < nqpts * ecmech::nsvec; i++) {
      stress_array[i] = 0.0;
   }

   //We'll leave these uninitialized for now since they're set in the
   //setup_data function.
   ddsdde_array = memoryManager::allocate<double>(nqpts * ecmech::nsvec * ecmech::nsvec);
   eng_int_array = memoryManager::allocate<double>(nqpts * ecmech::ne);
   w_vec_array = memoryManager::allocate<double>(nqpts * ecmech::nwvec);
   vol_ratio_array = memoryManager::allocate<double>(nqpts * ecmech::nvr);
   stress_svec_p_array = memoryManager::allocate<double>(nqpts * ecmech::nsvp);
   d_svec_p_array = memoryManager::allocate<double>(nqpts * ecmech::nsvp);

   std::cout << "Quat before: " << std::endl;
   for (int i = 0; i < ecmech::qdim; i++) {
      std::cout << state_vars[8 + i] << " ";
   }

   std::cout << std::endl;

   double stress_avg[6];
   double wts = 1.0 / nqpts;

   for (int i = 0; i < nsteps; i++) {
      std::cout << "About to setup data..." << std::endl;
      setup_data(device, nqpts, num_state_vars, dt, vgrad, stress_array, state_vars,
                 stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                 vol_ratio_array, eng_int_array);
      std::cout << "Data is now setup and now to run material model" << std::endl;
      //Now our kernel
      mat_model_kernel(device, mat_model_base, nqpts, dt,
                       num_state_vars, state_vars, stress_svec_p_array,
                       d_svec_p_array, w_vec_array, ddsdde_array,
                       vol_ratio_array, eng_int_array);
      std::cout << "Material model ran now to retrieve the data" << std::endl;
      //retrieve all of the data and put it back in the global arrays
      retrieve_data(device, nqpts, num_state_vars,
                    stress_svec_p_array, vol_ratio_array,
                    eng_int_array, state_vars, stress_array);
      for (int j = 0; j < ecmech::nsvec; j++) {
         stress_avg[j] = 0.0;
         for (int i_qpts = 0; i_qpts < nqpts; i_qpts++) {
            const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
            stress_avg[j] += wts * stress[j];
         }
      }

      std::cout << "Step# " << i + 1 << " Stress: ";
      for (int i = 0; i < ecmech::nsvec; i++) {
         std::cout << stress_avg[i] << " ";
      }

      std::cout << std::endl;
   }

   std::cout << "Quat final: ";
   for (int i = 0; i < ecmech::qdim; i++) {
      std::cout << state_vars[8 + i] << " ";
   }

   std::cout << std::endl;

   //for loop time
   //setup_data(...);
   // A stride-n index range [beg, end, n) using type int.
   // We might want to setup n here to be greater than one to potentially to have
   // multiple pts passed to our solver
   // int n = some number;
   // int num_threads = some number; When using CUDA or HIP we'll need to play
   // around with the number of threads used to reduce memory type pressures.
   // RAJA::TypedRangeStrideSegment<int> striden_range(0, nqpts, n);
   // model_kernel(...);
   // retrieve_data(...);

   //Delete all variables declared using new now.

   //Right now this is the only model that I'll be testing out.
   //After things work I'll add in the other models as well
#if defined(RAJA_ENABLE_CUDA)
   if (device == Accelerator::CUDA) {
      if (mat_model_str.compare("voce") == 0) {
         RAJA::forall<RAJA::cuda_exec<1> >(RAJA::RangeSegment(0, 1), [ = ] RAJA_HOST_DEVICE(int) {
            mat_modela->~matModel();
         });
      }
   }
#endif


   memoryManager::deallocate(mat_modela);
   memoryManager::deallocate(mat_modelb);
   memoryManager::deallocate(state_vars);
   memoryManager::deallocate(vgrad);
   memoryManager::deallocate(stress_array);
   memoryManager::deallocate(stress_svec_p_array);
   memoryManager::deallocate(d_svec_p_array);
   memoryManager::deallocate(w_vec_array);
   memoryManager::deallocate(ddsdde_array);
   memoryManager::deallocate(vol_ratio_array);
   memoryManager::deallocate(eng_int_array);
}
