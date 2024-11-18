#include "ECMech_cases.h"
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

   double dt = 0.00025;
   int nsteps = 60;

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
   std::vector<double> velocity_grad_init;
   //
   {
      // All the input arguments
      std::string option_file(argv[1]);

      std::string ori_file;
      std::string mat_prop_file;
      std::string device_type;
      std::string dt_vals = "0.00025";
      std::string nsteps_vals = "60";
      std::string velocity_grad_vals = "[[-0.5 0.0 0.0], [0.0 -0.5 0.0], [0.0 0.0 1.0]]";

      {
         std::ostringstream fail_str;
         fail_str << "Option file could not be correctly parsed." << std::endl
                  << "Option file contains: quat file path, material model, " << std::endl
                  << "material param file path, and device type each on their own line." << std::endl
                  << "Optionally, the option file past those required values can also contain:" << std::endl
                  << "dt value" << std::endl
                  << "number of steps value" << std::endl
                  << "velocity gradient as defined using the following notation [[# # #], [# # #], [# # #]]" << std::endl
                  << "Note each line in these optional values requires that the previous optional value also be defined"
                  << std::endl;

         std::ifstream ofile(option_file);
         ofile.clear();
         std::string line;

         std::getline(ofile, ori_file);
         std::getline(ofile, mat_model_str);
         std::getline(ofile, mat_prop_file);
         std::getline(ofile, device_type);

         if (ofile.fail()) {
            std::cerr << fail_str.str();
            return 1;
         }

         if (ofile.peek() != std::ifstream::traits_type::eof()) {
            std::getline(ofile, dt_vals);
         }
         if (ofile.peek() != std::ifstream::traits_type::eof()) {
            std::getline(ofile, nsteps_vals);
         }
         if (ofile.peek() != std::ifstream::traits_type::eof()) {
            std::getline(ofile, velocity_grad_vals);
         }

         if (dt_vals.size() == 0 || nsteps_vals.size() == 0 || velocity_grad_vals.size() == 0) {
            std::cerr << fail_str.str();
            std::cerr << "Check for an empty string for one of the optional variables" << std::endl;
            std::cerr << "dt_val.size()" << dt_vals.size()
                      << " nsteps_vals.size() " << nsteps_vals.size()
                      << " velocity_grad_vals.size() " << velocity_grad_vals.size() << std::endl;
            return 1;
         }
      }

      {
         std::istringstream iss(dt_vals);
         iss >> dt;
      }

      {
         std::istringstream iss(nsteps_vals);
         iss >> nsteps;
      }

      {
         std::istringstream iss(velocity_grad_vals);
         auto parse_data_row = [=] (auto& data, std::istringstream& stream) {
            constexpr auto max_size = std::numeric_limits<std::streamsize>::max();
            stream.ignore(max_size, '[');
            double vrow[3] = {};
            stream >> vrow[0] >> vrow[1] >> vrow[2];
            data.push_back(vrow[0]);
            data.push_back(vrow[1]);
            data.push_back(vrow[2]);
         };
         iss.ignore(1, '[');

         for (int i = 0; i < 3; i++) {
            parse_data_row(velocity_grad_init, iss);
         }
      }

      // This next chunk reads in all of the quaternions and pushes them to a vector.
      // It will exit if 4 values are not read on a line.
      bool quat_random = false;
      unsigned int quat_nrand = 1;
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
            // std::default_random_engine gen(42);
            // std::normal_distribution<double> distrib(0.0, 1.0); // An alternative way to initialize the quats
            // std::uniform_real_distribution<double> udistrib(-1.0, 1.0);
            std::minstd_rand0 gen(42);
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

      std::cout << "Orientation File: " << ori_file << std::endl;
      std::cout << "Material Property File: " << mat_prop_file << std::endl;
      std::cout << "Material Model: " << mat_model_str << std::endl;
      std::cout << "Execution Strategy: " << device_type << std::endl;
      std::cout << "Delta Time Step: " << dt << std::endl;
      std::cout << "Number of steps: " << nsteps << std::endl;
      std::cout << "Number of qpts: " << nqpts << std::endl;
      std::cout << "Velocity Gradient: " << std::endl;
      {
         auto it = velocity_grad_init.begin();
         for (int irow = 0; irow < 3; irow++) {
            for (int icol = 0; icol < 3; icol++) {
               std::cout << *it++ << " ";
            }
            std::cout << std::endl;
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

      std::cout << "num_props: " << num_props << " num_state_vars " << num_state_vars << std::endl;
      std::cout << "num_hardness: " << num_hardness << " num_gdot " << num_gdot << " iHistLbGdot " << iHistLbGdot << std::endl;

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
      // We're now initializing our state variables and velocity_grad to be used in other parts
      // of the simulations.
      constexpr size_t num_var_variables = (1 + ecmech::nsdd + + ecmech::ne + ecmech::nwvec + ecmech::nvr + ecmech::nsvec + 2 * ecmech::nsvp + ecmech::nsvec * ecmech::nsvec + ecmech::ndim * ecmech::ndim);
      const size_t num_items = nqpts * (num_state_vars + num_var_variables);
      auto mm = memoryManager<double>(num_items);
      auto state_vars = mm.getNew(nqpts * num_state_vars, class_device);
      auto velocity_grad = mm.getNew(nqpts * ecmech::ndim * ecmech::ndim, class_device);

      init_data(quats, mat_model_base, nqpts, num_hardness,
                num_gdot, iHistLbGdot, num_state_vars, state_vars);
      std::cout << "Data is now initialized" << std::endl;
      setup_velocity_grad(velocity_grad_init, velocity_grad, nqpts);

   // The stress array is the only one of the below variables that needs to be
   // initialized to 0.
   auto cauchy_stress_array = mm.getNew(nqpts * ecmech::nsvec, class_device);
   snls::forall(0, nqpts * ecmech::nsvec,
      [=]
      __ecmech_hdev__
      (int i) {
         cauchy_stress_array[i] = 0.0;
   });

   // We'll leave these uninitialized for now, since they're set in the
   // setup_data function.
   
   auto ddsdde_array = mm.getNew(nqpts * ecmech::nsvec * ecmech::nsvec, class_device);
   auto internal_energy_array = mm.getNew(nqpts * ecmech::ne, class_device);
   auto spin_vec_array = mm.getNew(nqpts * ecmech::nwvec, class_device);
   auto rel_vol_ratios_array = mm.getNew(nqpts * ecmech::nvr, class_device);
   auto cauchy_stress_d6p_array = mm.getNew(nqpts * ecmech::nsvp, class_device);
   auto def_rate_d6v_array = mm.getNew(nqpts * ecmech::nsvp, class_device);
   auto tkelv_array = mm.getNew(nqpts, class_device);
   auto sdd_array = mm.getNew(nqpts * ecmech::nsdd, class_device);

   double stress_avg[6];
   double wts = 1.0 / nqpts;

   RAJA::RangeSegment default_range(0, nqpts);

   RAJA::Timer run_time;

   run_time.start();

   for (int i = 0; i < nsteps; i++) {
      // set up our data in the correct format that the material model kernel expects
      setup_data(nqpts, num_state_vars, dt, velocity_grad, cauchy_stress_array, state_vars,
                 cauchy_stress_d6p_array, def_rate_d6v_array, spin_vec_array, ddsdde_array,
                 rel_vol_ratios_array, internal_energy_array, tkelv_array);
      // run our material model
      mat_model_kernel(mat_model_base, nqpts, dt,
                       state_vars, cauchy_stress_d6p_array,
                       def_rate_d6v_array, spin_vec_array, ddsdde_array,
                       rel_vol_ratios_array, internal_energy_array, tkelv_array, sdd_array);
      // retrieve all of the data and put it back in the global arrays
      retrieve_data(nqpts, num_state_vars,
                    cauchy_stress_d6p_array, rel_vol_ratios_array,
                    internal_energy_array, state_vars, cauchy_stress_array);

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
                  const double* cauchy_stress = &(cauchy_stress_array[i_qpts * ecmech::nsvec]);
                  seq_sum += wts * cauchy_stress[j];
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
                  const double* cauchy_stress = &(cauchy_stress_array[i_qpts * ecmech::nsvec]);
                  omp_sum += wts * cauchy_stress[j];
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
                  const double* cauchy_stress = &(cauchy_stress_array[i_qpts * ecmech::nsvec]);
                  gpu_sum += wts * cauchy_stress[j];
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
      // If we want to later output the deviatoric stress then we can add that in as
      // an option here with the following set of code...
      /*
      const double stress_mean = (stress_avg[0] + stress_avg[1] + stress_avg[2]) / 3.0;
      std::cout << "Deviatoric Stress: ";
      for (int i = 0; i < ecmech::ndim; i++) {
         std::cout << stress_avg[i] - stress_mean << " ";
      }
      for (int i = ecmech::ndim; i < ecmech::nsvec; i++) {
         std::cout << stress_avg[i] << " ";
      }
      std::cout << " " << stress_mean << std::endl;
      */
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

