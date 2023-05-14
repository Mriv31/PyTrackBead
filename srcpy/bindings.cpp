#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "tracking/Tracker.hpp"

void addBuffer(Tracker& vt, pybind11::array_t<unsigned char>& image, int image_size)
{
  pybind11::buffer_info buffer = image.request();
  vt.addBuffer((unsigned char*)buffer.ptr, image_size);
}

void define_image_counter(Tracker& vt, pybind11::array_t<int>& counter, int size)
{
  pybind11::buffer_info buffer = counter.request();
  vt.define_image_counter((int*)buffer.ptr, size);
}

void define_timebuffer(Tracker& vt, pybind11::array_t<uint64_t>& counter)
{
  pybind11::buffer_info buffer = counter.request();
  vt.define_timebuffer((uint64_t*)buffer.ptr);
}


MutexCV* get_mutex(Tracker& vt)
{
  MutexCV* mcv = vt.get_mutex();
  pybind11::set_shared_data("sharedmutex",mcv);
  return mcv;
}





PYBIND11_MODULE(bindings, m)
{

  m.doc() = "BeadTracker bindings";

  pybind11::class_<MutexCV,std::unique_ptr<MutexCV, pybind11::nodelete>>(m, "MutexCV",pybind11::module_local(false))
    .def(pybind11::init<>());


pybind11::class_<Tracker>(m, "Tracker")
  .def(pybind11::init<int, int>())
  .def("addBead", &Tracker::addBead)
  .def("addBuffer",&addBuffer)
  .def("define_image_counter",&define_image_counter)
  .def("Start",&Tracker::Start)
  .def("Stop",&Tracker::Stop)

  .def("get_x",&Tracker::get_x)
  .def("get_y",&Tracker::get_y)
  .def("set_saving",&Tracker::set_saving)
  .def("stop_saving",&Tracker::stop_saving)
  .def("define_timebuffer",&define_timebuffer)
  .def("get_copy_x_array",&Tracker::get_copy_x_array)
  .def("get_copy_y_array",&Tracker::get_copy_y_array)
  .def("set_frames_to_analyze",&Tracker::set_frames_to_analyze)
  .def("get_mutex",&Tracker::get_mutex)
  .def("notify_one",&Tracker::notify_one);



  }
