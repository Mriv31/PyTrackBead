
#include <thread>

#include "Tracker.hpp"
#include <cstdio>
#include <stdexcept>


Tracker::Tracker(int height,int width)
{
  ci = 0;
  w = width;
  h = height;
  nbead = 0;

  beadl = new Bead*[nbeadmax];

  buffer_list = (unsigned char**)calloc(MAX_TRACKING_BUFFER_SIZE,sizeof(unsigned char*));
  O_il = (O_i**)calloc(MAX_TRACKING_BUFFER_SIZE,sizeof(O_i*));

  buffer_size = 0;
  nb_frames_analyzed = 0;
  ci = 0;
  started = 0;
  saving = 0;
  defined_time_buffer = 0;
  stop = 0;
  stopped = 1;

  _frames_to_analyze = NULL;
  defined_image_counter = 0;

}

float Tracker::get_x(int bdnb)
{
  return beadl[bdnb]->get_x();

}
float Tracker::get_y(int bdnb)
{
  return beadl[bdnb]->get_y();

}




void Tracker::addBuffer(unsigned char* buffer, int size)
{
  if (size != w*h) throw std::runtime_error("Tried to add a Buffer with wrong size");
  if (started) throw std::runtime_error("Cannot add buffer if the tracker is started") ;
  if (buffer_size>=MAX_TRACKING_BUFFER_SIZE) throw std::runtime_error("Max buffer size reached"); // add runtime_error

  buffer_list[buffer_size] = buffer;
  O_il[buffer_size] = new O_i;
  O_il[buffer_size]->im.nx = w;
  O_il[buffer_size]->im.ny = h;
  O_il[buffer_size]->im.data_type = IS_CHAR_IMAGE;
  O_il[buffer_size]->im.pixel = new unsigned char*[h];
  for (int j=0;j<h;j++)
  {  O_il[buffer_size]->im.pixel[j] = &buffer[j * w];
     //printf("buf %d Oi[%d][5] = %d\n",buffer_size,j,O_il[buffer_size]->im.pixel[j][5]);
  }


  buffer_size+=1;
}

void Tracker::define_image_counter(int *buffer, int size)
{
  if (defined_image_counter) throw std::runtime_error("image counter already defined");
  if (size !=1) throw std::runtime_error("size of image counter must be 1");
  defined_image_counter = 1;

  _frames_to_analyze = buffer;
}

void Tracker::define_timebuffer(uint64_t *buffer)
{
  if (defined_time_buffer) throw std::runtime_error("image counter already defined");
  defined_time_buffer = 1;
  timebuffer = buffer;
}


void Tracker::addBead(int clh,int cwh, int xch, int ych)
{
  if (saving) throw std::runtime_error("cannot add bead during saving.");
  if (nbead >= nbeadmax) throw std::runtime_error("max number of beads exceeded");
  beadl[nbead] = new Bead(clh,cwh,xch,ych);
  nbead+=1;
}




void DoTrack(Tracker* tracker) {
  int frames_late = 0;
  int fta = 0;
  tracker->stopped = 0;
  tracker->stop = 0;

  std::unique_lock lock(tracker->mcv.m_mutex,std::defer_lock);
    while(1)
    {
      if (tracker->stop)
      {
        tracker->stop_saving();
        tracker->stopped = 1;
        return;
      }

      if (tracker->buffer_size == 0)
      {
        throw std::runtime_error("no buffer given to the track");
      }
      //here defined O_il
      lock.lock();
      tracker->mcv.cv.wait(lock);
      fta = *tracker->_frames_to_analyze;
      lock.unlock();

      frames_late = fta - tracker->nb_frames_analyzed;
      if (frames_late >= tracker->buffer_size)
      {
        throw std::runtime_error("missed more frames than buffer size. Will stop.");
      }
      while (frames_late > 0)
      {

        for (int i =0;i<tracker->nbead;i++)
        {

          tracker->beadl[i]->Track(tracker->O_il[tracker->ci]);

          if (tracker->saving)
          {

            tracker->outfile << tracker->timebuffer[tracker->ci] << ";" << tracker->beadl[i]->get_x() << ";" <<tracker->beadl[i]->get_y() ;
        }
          }

          if (tracker->saving)
          {

            tracker->outfile << "\n";
        }


      tracker->ci+=1;
      tracker->nb_frames_analyzed+=1;
      frames_late = fta - tracker->nb_frames_analyzed;

      if (tracker->ci == tracker->buffer_size) tracker->ci =0;
    }
     }

}

void Tracker::stop_saving()
{
  if (saving == 0) return ;
  outfile.close();
  saving = 0;
}



void Tracker::set_saving(std::string file)
{
  outfile = std::ofstream(file, std::ios::out);
  outfile << "Time";
  for (int i = 0;i<nbead;i++)
  {
    outfile << ";Bead" << i << "-x" << ";Bead" << i << "-y";
  }
  outfile << "\n";
  saving = 1;
}

void Tracker::Stop() {
    stop = 1;
    while(stopped == 0)
    {
      0;
    }
}
void Tracker::Start() {
    std::thread t(DoTrack, this);
    t.detach();
}

std::vector<float> Tracker::get_copy_x_array(int i)
{
  return beadl[i]->get_copy_x_array();
}


std::vector<float> Tracker::get_copy_y_array(int i)
{
  return beadl[i]->get_copy_y_array();
}
