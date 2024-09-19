#include "cases/ECMech_cases_fcc_defs.h"
#include "RAJA/RAJA.hpp"
#include "RAJA/util/Timer.hpp"
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
   // TODO:
   // Compare GPU versus the serial results (We need to figure out what the bounds on our
   // differences are between our results. We could see some due to CPU and GPU being different
   // devices with potentially subtle different ways of dealing with floating point arithmetic)

   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] <<
         " option file path which contains: quat file path, material model, "
                << "material param file path, and device type each on their own line."
                << std::endl;
      return 1;
   }

   const double dt = 0.00025;
   const int nsteps = 60;

   int nqpts = 0;
   int num_props = 0;
   int num_hardness = 0;
   int num_gdot = 0;
   int iHistLbGdot = 0;

   ecmech::matModelBase* mat_model_base;
   ecmech::ExecutionStrategy class_device;
   std::string mat_model_str;

   // The below scope of work sets up everything that we're going to be doing initially.
   // We're currently doing this to ensure that memory used during the set-up is freed
   // early on before we start doing all of our computations. We could probably keep everything
   // in scope without running into memory issues.
   //
   int num_state_vars;
   // Quaternion and the number of quaternions total.
   std::vector<double> quats;
   //
   {
      // All the input arguments
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

      // This next chunk reads in all of the quaternions and pushes them to a vector.
      // It will exit if 4 values are not read on a line.
      bool quat_random = false;
      unsigned int quat_nrand = 0;
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
            // provide a seed so things are reproducible
            std::default_random_engine gen(42);
            // std::normal_distribution<double> distrib(0.0, 1.0); // An alternative way to initialize the quats
            std::uniform_real_distribution<double> udistrib(-1.0, 1.0);
            std::vector<double> q_state = { 1., 0., 0., 0. };

            for (unsigned int i = 0; i < quat_nrand; i++) {
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
               } // error

               quats.push_back(q1); quats.push_back(q2); quats.push_back(q3); quats.push_back(q4);
            }
         }
      }


      // Read and store our material property data
      // We're going to check that the number of properties are what we expect
      // when we initialize the classes.
      std::vector<double> mat_props;
      int mp_nlines = 0;
      {
         std::ifstream mfile(mat_prop_file);
         std::string line;
         while (std::getline(mfile, line)) {
            std::istringstream iss(line);
            double prop;

            mp_nlines += 1;
            if (!(iss >> prop)) {
               std::cerr << "Material prop file has a malformed line on line: " << mp_nlines << std::endl;
               return 1;
            } // error

            mat_props.push_back(prop);
         }
      }


      // We now detect which device is desired to run the different cases.
      // Compiler flags passed in will tell us which options are available
      // based on what RAJA was built with. If we do not have support for
      // the chosen value then we should error out and let the user know
      // which values are available.

      if (device_type.compare("CPU") == 0) {
         class_device = ECM_EXEC_STRAT_CPU;
      }
#if defined(RAJA_ENABLE_OPENMP)
      else if (device_type.compare("OpenMP") == 0) {
         class_device = ECM_EXEC_STRAT_OPENMP;
      }
#endif
#if defined(RAJA_ENABLE_CUDA) || defined(RAJA_ENABLE_HIP)
      else if (device_type.compare("GPU") == 0) {
         class_device = ECM_EXEC_STRAT_GPU;
      }
#endif
      else {
         std::cerr << "Accelerator is not supported or RAJA was not built with" << std::endl;
         return 1;
      }

      // Some basic set up for when we initialize our class
      // Opts and strs are just empty vectors of int and strings
      std::vector<double> params;
      std::vector<int> opts;
      std::vector<std::string> strs;

      for (unsigned int i = 0; i < mat_props.size(); i++) {
         params.push_back(mat_props.at(i));
      }

      std::cout << "\nAbout to initialize class" << std::endl;
      mat_model_base = ecmech::makeMatModel(mat_model_str);
      auto index_map = ecmech::modelParamIndexMap(mat_model_str);
      num_props = index_map["num_params"];
      num_state_vars = index_map["num_hist"];
      num_state_vars += ecmech::ne + 1;

      num_hardness = index_map["num_hardening"];
      num_gdot = index_map["num_slip_system"];
      iHistLbGdot = index_map["index_slip_rates"];

      // This check used to be in the loop used to read in the material properties
      // However, things were re-arranged, so it's now during the class initialization
      if (mp_nlines != num_props) {
         std::cerr << "Material prop file should have " << num_props
                     << " properties (each on their own line). A total of " << mp_nlines
                     << " properties were provided instead." << std::endl;
         return 1;
      }

      std::vector<size_t> strides;
      // Deformation rate stride
      strides.push_back(ecmech::nsvp);
      // Spin rate stride
      strides.push_back(ecmech::ndim);
      // Volume ratio stride
      strides.push_back(ecmech::nvr);
      // Internal energy stride
      strides.push_back(ecmech::ne);
      // Stress vector stride
      strides.push_back(ecmech::nsvp);
      // History variable stride
      strides.push_back(num_state_vars);
      // Temperature stride
      strides.push_back(1);
      // SDD stride
      strides.push_back(ecmech::nsdd);

      mat_model_base->updateStrides(strides);

      // We really shouldn't see this change over time at least for our applications.
      mat_model_base->setExecutionStrategy(class_device);
      mat_model_base->initFromParams(opts, params, strs);
      mat_model_base->complete();

      std::cout << "Class has been completely initialized" << std::endl;
   }
      // We're now initializing our state variables and vgrad to be used in other parts
      // of the simulations.
      constexpr size_t num_var_variables = (1 + ecmech::nsdd + + ecmech::ne + ecmech::nwvec + ecmech::nvr + ecmech::nsvec + 2 * ecmech::nsvp + ecmech::nsvec * ecmech::nsvec + ecmech::ndim * ecmech::ndim);
      const size_t num_items = nqpts * (num_state_vars + num_var_variables);
      auto mm = memoryManager<double>(num_items);
      auto state_vars = mm.getNew(nqpts * num_state_vars, class_device);
      auto vgrad = mm.getNew(nqpts * ecmech::ndim * ecmech::ndim, class_device);

      init_data(quats, mat_model_base, nqpts, num_hardness,
                num_gdot, iHistLbGdot, num_state_vars, state_vars);
      std::cout << "Data is now initialized" << std::endl;
      setup_vgrad(vgrad, nqpts);

   // The stress array is the only one of the below variables that needs to be
   // initialized to 0.
   auto stress_array = mm.getNew(nqpts * ecmech::nsvec, class_device);
   snls::forall(0, nqpts * ecmech::nsvec,
      [=]
      __ecmech_hdev__
      (int i) {
         stress_array[i] = 0.0;
   });

   // We'll leave these uninitialized for now, since they're set in the
   // setup_data function.
   
   auto ddsdde_array = mm.getNew(nqpts * ecmech::nsvec * ecmech::nsvec, class_device);
   auto eng_int_array = mm.getNew(nqpts * ecmech::ne, class_device);
   auto w_vec_array = mm.getNew(nqpts * ecmech::nwvec, class_device);
   auto vol_ratio_array = mm.getNew(nqpts * ecmech::nvr, class_device);
   auto stress_svec_p_array = mm.getNew(nqpts * ecmech::nsvp, class_device);
   auto d_svec_p_array = mm.getNew(nqpts * ecmech::nsvp, class_device);
   auto temp_array = mm.getNew(nqpts, class_device);
   auto sdd_array = mm.getNew(nqpts * ecmech::nsdd, class_device);

   double stress_avg[6];
   double wts = 1.0 / nqpts;

   RAJA::RangeSegment default_range(0, nqpts);

   RAJA::Timer run_time;

   run_time.start();

   for (int i = 0; i < nsteps; i++) {
      // set up our data in the correct format that the material model kernel expects
      setup_data(nqpts, num_state_vars, dt, vgrad, stress_array, state_vars,
                 stress_svec_p_array, d_svec_p_array, w_vec_array, ddsdde_array,
                 vol_ratio_array, eng_int_array, temp_array);
      // run our material model
      mat_model_kernel(mat_model_base, nqpts, dt,
                       state_vars, stress_svec_p_array,
                       d_svec_p_array, w_vec_array, ddsdde_array,
                       vol_ratio_array, eng_int_array, temp_array, sdd_array);
      // retrieve all of the data and put it back in the global arrays
      retrieve_data(nqpts, num_state_vars,
                    stress_svec_p_array, vol_ratio_array,
                    eng_int_array, state_vars, stress_array);

      switch ( class_device ) {
         default :
         case ECM_EXEC_STRAT_CPU :
         {
            if (NEVALS_COUNTS) {
               RAJA::ReduceSum<RAJA::seq_reduce, double> seq_sum(0.0);
               RAJA::ReduceMin<RAJA::seq_reduce, double> seq_min(100.0); // We know this shouldn't ever be more than 100
               RAJA::ReduceMax<RAJA::seq_reduce, double> seq_max(0.0); // We know this will always be at least 1.0
               RAJA::forall<RAJA::seq_exec>(default_range, [ = ] (int i_qpts){
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
               RAJA::forall<RAJA::seq_exec>(default_range, [ = ] (int i_qpts){
                  const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
                  seq_sum += wts * stress[j];
               });
               stress_avg[j] = seq_sum.get();
	    }
         }
         break;
#if defined(RAJA_ENABLE_OPENMP)
         case ECM_EXEC_STRAT_OPENMP :
         {   
            if (NEVALS_COUNTS) {
               RAJA::ReduceSum<RAJA::omp_reduce_ordered, double> omp_sum(0.0);
               RAJA::ReduceMin<RAJA::omp_reduce_ordered, double> omp_min(100.0); // We know this shouldn't ever be more than 100
               RAJA::ReduceMax<RAJA::omp_reduce_ordered, double> omp_max(0.0); // We know this will always be at least 1.0
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
         break;
#endif
#if defined(RAJA_ENABLE_CUDA) || defined(RAJA_ENABLE_HIP)
         case ECM_EXEC_STRAT_GPU :
         {
#if defined(RAJA_ENABLE_CUDA)
            using gpu_reduce = RAJA::cuda_reduce;
            using gpu_policy = RAJA::cuda_exec<1024>;
#else
            using gpu_reduce = RAJA::hip_reduce;
            using gpu_policy = RAJA::hip_exec<1024>;
#endif
            if (NEVALS_COUNTS) {
               RAJA::ReduceSum<gpu_reduce, double> gpu_sum(0.0);
               RAJA::ReduceMin<gpu_reduce, double> gpu_min(100.0); // We know this shouldn't ever be more than 100
               RAJA::ReduceMax<gpu_reduce, double> gpu_max(0.0); // We know this will always be at least 1.0
               RAJA::forall<gpu_policy>(default_range, [ = ] RAJA_DEVICE(int i_qpts){
                  double* nfunceval = &(state_vars[i_qpts * num_state_vars + 2]);
                  gpu_sum += wts * nfunceval[0];
                  gpu_max.max(nfunceval[0]);
                  gpu_min.min(nfunceval[0]);
               });
               std::cout << "Min Func Eval: " << gpu_min.get() << " Mean Func Evals: " <<
                  gpu_sum.get() << " Max Func Eval: " << gpu_max.get() << std::endl;
            }
            for (int j = 0; j < ecmech::nsvec; j++) {
               RAJA::ReduceSum<gpu_reduce, double> gpu_sum(0.0);
               RAJA::forall<gpu_policy>(default_range, [ = ] RAJA_DEVICE(int i_qpts){
                  const double* stress = &(stress_array[i_qpts * ecmech::nsvec]);
                  gpu_sum += wts * stress[j];
               });
               stress_avg[j] = gpu_sum.get();
            }
         }
         break;
#endif
      } // switch ( class_device ) 

      // On CORAL architectures these print statements don't really add anything to the execution time.
      // So, we're going to keep them to make sure things are correct between the different runs.
      std::cout << "Step# " << i + 1 << " Stress: ";
      for (int i = 0; i < ecmech::nsvec; i++) {
         std::cout << stress_avg[i] << " ";
      }

      std::cout << std::endl;
   }

   run_time.stop();

   double time = run_time.elapsed();

   std::cout << std::endl;

   std::cout << "Run time of set-up, material, and retrieve kernels over " <<
      nsteps << " time steps is: " << time << "(s)" << std::endl;
   // All the variables share the same memory buffer so once the mm object goes out of scope
   // it's deconstructor will free all of the memory used

   return 0;
}

