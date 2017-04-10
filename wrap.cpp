//
//  wrapper.cpp
//  PhD
//
//  Created by Louis Ponet on 25/12/2016.
//  Copyright © 2016 Louis Ponet. All rights reserved.
//

//#include "calc.hpp"
//#include "ham_calc.hpp"
//#include "sup_calc.hpp"
//#include "dip_calc.hpp"
#include "typedefs.hpp"

namespace py = pybind11;

void init_calc(py::module& m);
void init_sup_calc(py::module& m);
void init_dip_calc(py::module& m);
void init_ham_calc(py::module& m);

PYBIND11_PLUGIN(wrap) {
    py::module m("wrap", "pybind11 example plugin");
    init_calc(m);
    init_sup_calc(m) ;
    init_ham_calc(m);
    
    init_dip_calc(m);
    py::enum_<AngMomentum>(m,"AngMomentum").value("Lx",AngMomentum::Lx).value("Ly",AngMomentum::Ly).value("Lz",AngMomentum::Lz).export_values();
    py::enum_<Direction>(m,"Direction").value("Ex",Direction::Ex).value("Ey",Direction::Ey).value("Ez",Direction::Ez).value("Ex_",Direction::Ex_).value("Ey_",Direction::Ey_).value("Ez_",Direction::Ez_).value("END",Direction::END).export_values();
    return m.ptr();
}
