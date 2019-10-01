#include "ECMech_cases.h"
#include "ECMech_evptnWrap.h"
#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
#if defined(RAJA_ENABLE_CUDA)
#include "cuda.h"
//#include "cuda_profiler_api.h"
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

#define NEVALS_COUNTS false


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

   std::vector<uint> strides;
   strides.push_back(ecmech::nsvp);
   strides.push_back(ecmech::ndim);
   strides.push_back(ecmech::nvr);
   strides.push_back(ecmech::ne);
   strides.push_back(ecmech::nsvp);
   strides.push_back(num_state_vars);
   strides.push_back(1);
   strides.push_back(ecmech::nsdd);

   ecmech::matModelEvptn_FCC_A mat_modela(strides.data(), strides.size());
   ecmech::matModelEvptn_FCC_B mat_modelb(strides.data(), strides.size());

   ecmech::Accelerator class_device;

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
   double* temp_array = nullptr;
   double* sdd_array = nullptr;

   //If this is turned off then we do take a performance hit. At least if this is 10k and above things
   //run fine.
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
      bool quat_random = false;
      uint quat_nrand = 0;
      {
         std::ifstream qfile(ori_file);
         std::string line;
         {
            std::getline(qfile, line);
            std::istringstream iss(line);
            std::string tmp_str;

            if (!(iss >> tmp_str >> quat_nrand)) {
               std::cerr << "Quat file starting line should be either the following in parantheses: " <<
                  "(#random num_quats) where num_quats is a positive value for the number of " <<
                  "quaternions that you want randomly generated, or it can be (#data 0) where " <<
                  "reads in all of the content of the file where each line is a quat" << std::endl;
               return 1;
            }

            if (tmp_str.compare("#random") == 0) {
               quat_random = true;
               nqpts = quat_nrand;
            }
         }
         if (quat_random) {
            //provide a seed so things are reproducible
            std::default_random_engine gen(42);
            // std::normal_distribution<double> distrib(0.0, 1.0);
            std::uniform_real_distribution<double> udistrib(-1.0, 1.0);
            std::vector<double> q_state = { 1., 0., 0., 0. };

            for (int i = 0; i < quat_nrand; i++) {
               q_state[0] = udistrib(gen);
               q_state[1] = udistrib(gen);
               q_state[2] = udistrib(gen);
               q_state[3] = udistrib(gen);

               ecmech::vecsVNormalize<ecmech::qdim>(q_state.data());

               quats.push_back(q_state[0]);
               quats.push_back(q_state[1]);
               quats.push_back(q_state[2]);
               quats.push_back(q_state[3]);
            }
         }
         else {
            while (std::getline(qfile, line)) {
               std::istringstream iss(line);
               double q1, q2, q3, q4;

               nqpts += 1;

               if (!(iss >> q1 >> q2 >> q3 >> q4)) {
                  std::cerr << "Quat file has malformed line on line: " << nqpts << std::endl;
                  return 1;
               }// error

               quats.push_back(q1); quats.push_back(q2); quats.push_back(q3); quats.push_back(q4);
            }
         }
      }


      //Read and store our material property data
      //We're going to check that the number of properties are what we expect
      //when we initialize the classes.
      std::vector<double> mat_props;
      int mp_nlines = 0;
      {
         std::ifstream mfile(mat_prop_file);
         std::string line;
         while (std::getline(mfile, line)) {
            std::istringstream iss(line);
            double prop;

            mp_nlines += 1;
            // std::cout << line << std::endl;
            if (!(iss >> prop)) {
               std::cerr << "Material prop file has a malformed line on line: " << mp_nlines << std::endl;
               return 1;
            }    // error

            mat_props.push_back(prop);
         }
      }


      //We now detect which device is desired to run the different cases.
      //Compiler flags passed in will tell us which options are available
      //based on what RAJA was built with. If we do not have support for
      //the chosen value then we should error out and let the user know
      //which values are available.

      if (device_type.compare("CPU") == 0) {
         class_device = ecmech::Accelerator::CPU;
      }
      #if defined(RAJA_ENABLE_OPENMP)
      else if (device_type.compare("OpenMP") == 0) {
         class_device = ecmech::Accelerator::OPENMP;
      }
      #endif
      #if defined(RAJA_ENABLE_CUDA)
      else if (device_type.compare("CUDA") == 0) {
         class_device = ecmech::Accelerator::CUDA;
      }
      #endif
      else {
         std::cerr << "Accelerator is not supported or RAJA was not built with" << std::endl;
         return 1;
      }

      //Some basic set up for when we initialize our class
      //Opts and strs are just empty vectors of int and strings
      std::vector<double> params;
      std::vector<int> opts;
      std::vector<std::string> strs;

      for (uint i = 0; i < mat_props.size(); i++) {
         params.push_back(mat_props.at(i));
      }

      std::cout << "\nAbout to initialize class" << std::endl;
      //Initialize our base class using the appropriate model
      if (mat_model_str.compare("voce") == 0) {
         num_props = 17;
         num_hardness = mat_modela.nH;
         num_gdot = mat_modela.nslip;
         iHistLbGdot = mat_modela.iHistLbGdot;

         //This check used to be in the loop used to read in the material properties
         //However, if things we're re-arranged, so it's now during the class initialization
         if (mp_nlines != num_props) {
            std::cerr << "Material prop file should have " << num_props
                      << " properties (each on their own line). A total of " << mp_nlines
                      << " properties were provided instead." << std::endl;
            return 1;
         }

         //We really shouldn't see this change over time at least for our applications.
         mat_modela.initFromParams(opts, params, strs, class_device);
         mat_modela.complete();
         mat_model_base = dynamic_cast<matModelBase*>(&mat_modela);
      }
      else if (mat_model_str.compare("mts") == 0) {
         num_props = 24;
         num_hardness = mat_modelb.nH;
         num_gdot = mat_modelb.nslip;
         iHistLbGdot = mat_modelb.iHistLbGdot;

         //This check used to be in the loop used to read in the material properties
         //However, if things we're re-arranged, so it's now during the class initialization
         if (mp_nlines != num_props) {
            std::cerr << "Material prop file should have " << num_props
                      << " properties (each on their own line). A total of " << mp_nlines
                      << " properties were provided instead." << std::endl;
            return 1;
         }

         //We really shouldn't see this change over time at least for our applications.
         mat_modelb.initFromParams(opts, params, strs, class_device);
         //
         mat_modelb.complete();
         mat_model_base = dynamic_cast<matModelBase*>(&mat_modelb);
      }
      else {
         std::cerr << "material model must be either voce or mts " << std::endl;
         return 1;
      }

      std::cout << "Class has been completely initialized" << std::endl;

      //We're now initializing our state variables and vgrad to be used in other parts
      //of the simulations.
      state_vars = memoryManager::allocate<double>(num_state_vars * nqpts);
      vgrad = memoryManager::allocate<double>(nqpts * ecmech::ndim * ecmech::ndim);

      double* quats_array = quats.data();

      init_data(class_device, quats_array, mat_model_base, nqpts, num_hardness,
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
   temp_array = memoryManager::allocate<double>(nqpts);
   sdd_array = memoryManager::allocate<double>(nqpts * ecmech::nsdd);

   std::cout << std::endl;

   double stress_avg[6];
   double wts = 1.0 / nqpts;

   RAJA::RangeSegment default_range(0, nqpts);

   RAJA::Timer run_time;

   run_time.start();

   // For profiling uses
   // #if defined(RAJA_ENABLE_CUDA)
   //    cudaProfilerStart();
   // #endif

   for (int i = 0; i < nsteps; i++) {
      // std::cout << "About to setup data..." << std::endl;
      setup_data(class_device, nqpts, num_state_vars, dt, vgrad, stress_array, state_vars,
                 stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                 vol_ratio_array, eng_int_array, temp_array);
      // std::cout << "Data is now setup and now to run material model" << std::endl;
      //Now our kernel
      mat_model_kernel(mat_model_base, nqpts, dt,
                       num_state_vars, state_vars, stress_svec_p_array,
                       d_svec_p_array, w_vec_array, ddsdde_array,
                       vol_ratio_array, eng_int_array, temp_array, sdd_array);
      // std::cout << "Material model ran now to retrieve the data" << std::endl;
      //retrieve all of the data and put it back in the global arrays
      retrieve_data(class_device, nqpts, num_state_vars,
                    stress_svec_p_array, vol_ratio_array,
                    eng_int_array, state_vars, stress_array);

      if (class_device == ecmech::Accelerator::CPU) {
         if (NEVALS_COUNTS) {
            RAJA::ReduceSum<RAJA::seq_reduce, double> seq_sum(0.0);
            RAJA::ReduceMin<RAJA::seq_reduce, double> seq_min(100.0);//We now this shouldn't ever be more than 100
            RAJA::ReduceMax<RAJA::seq_reduce, double> seq_max(0.0);//We now this will always be at least 1.0
            RAJA::forall<RAJA::loop_exec>(default_range, [ = ] (int i_qpts){
               double* nfunceval = &(state_vars[i_qpts * num_state_vars + 2]);
               seq_sum += wts * nfunceval[0];
               seq_max.max(nfunceval[0]);
               seq_min.min(nfunceval[0]);
            });
            std::cout << "Min Func Eval: " << seq_min.get() << " Mean Func Evals: " <<
               seq_sum.get() << " Max Func Eval: " << seq_max.get() << std::endl;
         }
         for (int j = 0; j < ecmech::nsvec; j++) {
            RAJA::ReduceSum<RAJA::seq_reduce, double> seq_sum(0.0);
            RAJA::forall<RAJA::loop_exec>(default_range, [ = ] (int i_qpts){
               const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
               seq_sum += wts * stress[j];
            });
            stress_avg[j] = seq_sum.get();
         }
      }
      #if defined(RAJA_ENABLE_OPENMP)
      if (class_device == ecmech::Accelerator::OPENMP) {
         if (NEVALS_COUNTS) {
            RAJA::ReduceSum<RAJA::omp_reduce_ordered, double> omp_sum(0.0);
            RAJA::ReduceMin<RAJA::omp_reduce_ordered, double> omp_min(100.0);//We now this shouldn't ever be more than 100
            RAJA::ReduceMax<RAJA::omp_reduce_ordered, double> omp_max(0.0);//We now this will always be at least 1.0
            RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i_qpts){
               double* nfunceval = &(state_vars[i_qpts * num_state_vars + 2]);
               omp_sum += wts * nfunceval[0];
               omp_max.max(nfunceval[0]);
               omp_min.min(nfunceval[0]);
            });
            std::cout << "Min Func Eval: " << omp_min.get() << " Mean Func Evals: " <<
               omp_sum.get() << " Max Func Eval: " << omp_max.get() << std::endl;
         }
         for (int j = 0; j < ecmech::nsvec; j++) {
            RAJA::ReduceSum<RAJA::omp_reduce_ordered, double> omp_sum(0.0);
            RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i_qpts){
               const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
               omp_sum += wts * stress[j];
            });
            stress_avg[j] = omp_sum.get();
         }
      }
      #endif
      #if defined(RAJA_ENABLE_CUDA)
      if (class_device == ecmech::Accelerator::CUDA) {
         if (NEVALS_COUNTS) {
            RAJA::ReduceSum<RAJA::cuda_reduce, double> cuda_sum(0.0);
            RAJA::ReduceMin<RAJA::cuda_reduce, double> cuda_min(100.0);//We now this shouldn't ever be more than 100
            RAJA::ReduceMax<RAJA::cuda_reduce, double> cuda_max(0.0);//We now this will always be at least 1.0
            RAJA::forall<RAJA::cuda_exec<1024> >(default_range, [ = ] RAJA_DEVICE(int i_qpts){
               double* nfunceval = &(state_vars[i_qpts * num_state_vars + 2]);
               cuda_sum += wts * nfunceval[0];
               cuda_max.max(nfunceval[0]);
               cuda_min.min(nfunceval[0]);
            });
            std::cout << "Min Func Eval: " << cuda_min.get() << " Mean Func Evals: " <<
               cuda_sum.get() << " Max Func Eval: " << cuda_max.get() << std::endl;
         }
         for (int j = 0; j < ecmech::nsvec; j++) {
            RAJA::ReduceSum<RAJA::cuda_reduce, double> cuda_sum(0.0);
            RAJA::forall<RAJA::cuda_exec<1024> >(default_range, [ = ] RAJA_DEVICE(int i_qpts){
               const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
               cuda_sum += wts * stress[j];
            });
            stress_avg[j] = cuda_sum.get();
         }
      }
      #endif

      //On Summit like architectures these print statements don't really add anything to the execution time.
      //So, we're going to keep them to make things are correct.
      std::cout << "Step# " << i + 1 << " Stress: ";
      for (int i = 0; i < ecmech::nsvec; i++) {
         std::cout << stress_avg[i] << " ";
      }

      std::cout << std::endl;
   }

   // For profiling uses
   // #if defined(RAJA_ENABLE_CUDA)
   //    cudaProfilerStart();
   // #endif

   run_time.stop();

   double time = run_time.elapsed();

   std::cout << std::endl;

   std::cout << "Run time of set-up, material, and retrieve kernels over " <<
      nsteps << " time steps is: " << time << "(s)" << std::endl;

   //Delete all variables declared using new now.

   memoryManager::deallocate(state_vars);
   memoryManager::deallocate(vgrad);
   memoryManager::deallocate(stress_array);
   memoryManager::deallocate(stress_svec_p_array);
   memoryManager::deallocate(d_svec_p_array);
   memoryManager::deallocate(w_vec_array);
   memoryManager::deallocate(ddsdde_array);
   memoryManager::deallocate(vol_ratio_array);
   memoryManager::deallocate(eng_int_array);
   memoryManager::deallocate(temp_array);
   memoryManager::deallocate(sdd_array);

   return 0;
}
