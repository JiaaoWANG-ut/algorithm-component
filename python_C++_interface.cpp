python3 链接C++ 使用pybind11,反复调用
======================================

//pip install pybind11
  
======================================example.cpp
#include <pybind11/pybind11.h>

int add(int a, int b) {
    return a + b;
}

namespace py = pybind11;

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin";
    m.def("add", &add, "A function which adds two numbers");
}
======================================

//pip install pybind11
//执行下面的命令
c++ -O3 -Wall -shared -std=c++11 -fPIC `python3-config --includes` example.cpp -o example`python3-config --extension-suffix`

======================================1.py

import example

# 反复调用 C++ 函数
for i in range(10000):
    result = example.add(3, 4)
    print(f"3 + 4 = {result}")
      
======================================
      
python3 1.py
