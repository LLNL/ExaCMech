#define DUMPVEC(aname, a) std::cout << "# " << aname << " : "; for (unsigned int iThing = 0; \
                                                                    iThing<a.size(); \
                                                                           ++iThing) { if (iThing) { std::cout << ","; \
                                                                                       } std::cout << a[iThing]; \
   } std::cout << std::endl;

double density0 = 3.0, cvav = 2.0e-5;
double tolerance = 1e-10;

