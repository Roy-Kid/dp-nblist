
#include "cuda_nblist.h"
#include "box.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(cuda_nblist_py, m) {
    py::class_<dpnblist::Box>(m, "Box")
        .def(py::init<std::vector<float>, std::vector<float>>(), py::arg("lengths"), py::arg("angles") = std::make_tuple(90.0f, 90.0f, 90.0f));

    py::class_<dpnblist::base_NBL>(m, "base_NBL")
        .def("build", &dpnblist::base_NBL::build)
        .def("update", &dpnblist::base_NBL::update)
        .def("gettime", &dpnblist::base_NBL::gettime)
        .def("out", &dpnblist::base_NBL::out)
        .def("get_listArray", &dpnblist::base_NBL::get_listArray);

    py::class_<dpnblist::CudaCellList>(m, "CellList")
        //.def(py::init<std::vector<float>&, float, float>())
        .def(py::init<dpnblist::Box*, float, float>())
        .def("build", &dpnblist::CudaCellList::build)
        .def("update", &dpnblist::CudaCellList::update)
        .def("gettime", &dpnblist::CudaCellList::gettime)
        .def("out", &dpnblist::CudaCellList::out)
        .def("get_listArray", &dpnblist::CudaCellList::get_listArray);

    
    m.def("createNeighborList", &dpnblist::NeighborList::createNeighborList, "Create a NeighborList object");
}
