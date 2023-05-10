#pragma once
#include "Bead.hpp"
# include "track_util.h"
#include <iostream>
#include <fstream>

#define MAX_TRACKING_BUFFER_SIZE 1024



class Tracker {
public:
  Tracker(int width,int length);
  void addBead(int clh,int cwh, int xch, int ych);
  void startBeadTrack();
  void addBuffer(unsigned char* buffer, int size);
  void Start();
  void Stop();
  void define_image_counter(long int *buffer,int size);
  float get_x(int bdnb);
  float get_y(int bdnb);
  void set_saving(std::string file);
  void define_timebuffer(uint64_t *buffer);
  std::vector<float> get_copy_x_array(int i);
  std::vector<float> get_copy_y_array(int i);
  void set_frames_to_analyze(long int a){*_frames_to_analyze=a;};





  int ci;
  long int *_frames_to_analyze;
  long int nb_frames_analyzed;
  int buffer_size;
  int beadtracking;
  O_i ** O_il;
  int nbead;
  Bead **beadl;
  unsigned char ** buffer_list;
  std::ofstream outfile;
  int saving;
  uint64_t* timebuffer;
  int defined_time_buffer;
  int stop;



private:
  int w;
  int h;
  unsigned char **buffer;


  const int nbeadmax = 500;
  int started;
  int defined_image_counter;



};
