#include "lbm.hpp"

int main() {
    LBM lbm(D2Q9, SRT, VELOCITY); 
    lbm.run();
    return 0;
}
