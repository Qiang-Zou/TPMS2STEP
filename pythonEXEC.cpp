#include "pythonEXEC.h"

void pythonEXEC(std::string modelType, double offsetValue1, double offsetValue2) {
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append('./')");
    std::string input = "python3 ../python/generateSFP_" + modelType + ".py -f" + to_string(offsetValue1);
    system(input.c_str());
    input = "python3 ../python/generateSFP_" + modelType + ".py -f" + to_string(offsetValue1 * (-1));
    system(input.c_str());
    if (offsetValue1 + offsetValue2 != 0) {
        input= "python3 ../python/generateSFP_" + modelType + ".py -f" + to_string(offsetValue2);
        system(input.c_str());
        input = "python3 ../python/generateSFP_" + modelType + ".py -f" + to_string(offsetValue2 * (-1));
        system(input.c_str());
    }
    Py_Finalize();
    std::cout << "TPMS2STEP > Generate SFP (Python) ok!" << endl;
}
