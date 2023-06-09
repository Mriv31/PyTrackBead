#pragma once
#include "Bead.hpp"
# include "track_util.h"
#include <iostream>
#include <fstream>
#include <mutex>
#include <condition_variable>

#define MAX_TRACKING_BUFFER_SIZE 1024


class MutexCV {
public:
    MutexCV() = default;
    void notify_one() { cv.notify_one(); };
    std::mutex m_mutex;
    std::condition_variable cv;


};


class Tracker {
public:
  Tracker(int width,int length);
  void addBead(int clh,int cwh, int xch, int ych);
  void addBuffer(unsigned char* buffer, int size);
  void Start();
  void Stop();
  void define_image_counter(int *buffer,int size);
  float get_x(int bdnb);
  float get_y(int bdnb);
  void set_saving(std::string file);
  void stop_saving();
  void define_timebuffer(uint64_t *buffer);
  std::vector<float> get_copy_x_array(int i);
  std::vector<float> get_copy_y_array(int i);
  void set_frames_to_analyze(int a){*_frames_to_analyze=a;};
  void notify_one() {mcv.notify_one();};
  MutexCV* get_mutex(){return &mcv;};





  int ci;
  int *_frames_to_analyze;
  int nb_frames_analyzed;
  int buffer_size;
  O_i ** O_il;
  int nbead;
  Bead **beadl;
  unsigned char ** buffer_list;
  std::ofstream outfile;
  int saving;
  uint64_t* timebuffer;
  int defined_time_buffer;
  int stop;
  int stopped;

  MutexCV mcv;



private:
  int w;
  int h;
  unsigned char **buffer;


  const int nbeadmax = 500;
  int started;
  int defined_image_counter;



};
