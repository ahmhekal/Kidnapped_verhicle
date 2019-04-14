/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */
#include "particle_filter.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "helper_functions.h"
using namespace std;
#define pi 3.14159265359

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {


  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
    num_particles = 100;  // TODO: Set the number of particles

    normal_distribution<double> dist_x (0, std[0]);
    normal_distribution<double> dist_y (0, std[1]);
    normal_distribution<double> dist_theta (0, std[2]);

  /*
    normal_distribution<double> noise_x(0, std[0]);
  	normal_distribution<double> noise_y(0, std[1]);
  	normal_distribution<double> noise_theta(0, std[2]);
  */
  Particle p;
    for (unsigned int i = 0; i < num_particles; ++i) {
    		double sample_x, sample_y, sample_theta;
        sample_x = dist_x(gen);
      	sample_y = dist_y(gen);
      	sample_theta = dist_theta(gen);
			  p.id=i;
      	p.x=x+sample_x;
        p.y=y+sample_y;
        p.theta=theta+sample_theta;

  /*    //Adding noise
      		p.x+=noise_x;
      		p.y+=noise_y;
      		p.theta+=noise_theta;
 */
			  p.weight=1.0;
      	particles.push_back(p);
    }
    is_initialized = true;

}


void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
                                  default_random_engine gen;

    /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  double new_x, new_y,new_theta,theta;

  if (fabs(yaw_rate)>0.0001){
     for (unsigned int i=0;i<num_particles;++i)
  		{
        theta=particles[i].theta;
        new_theta=yaw_rate*delta_t;

        new_x=(velocity/yaw_rate)*( sin(theta+new_theta) - sin(theta) );
        new_y=(velocity/yaw_rate)*( cos(theta) - cos(theta+new_theta) );

        particles[i].x+=new_x;
        particles[i].y+=new_y;
        particles[i].theta+=new_theta;
     }
  }

  else
  {
    for (unsigned int i=0;i<num_particles;++i)
  		{
       particles[i].x += velocity * delta_t * cos(particles[i].theta);
       particles[i].y += velocity * delta_t * sin(particles[i].theta);
      }
  }

  //Adding noise
      for (unsigned int i=0;i<num_particles;++i)
  		{

        normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);

      }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */


    for (unsigned int i=0; i<observations.size(); ++i)
    {
      int nearest=0;
      double mindist=numeric_limits<float>::max();

      for(unsigned int j=0;j<predicted.size();++j)
      {
        double distance = dist(predicted[j].x,predicted[j].y,observations[i].x,observations[i].y);
        if(distance<mindist) //if it is nearest
          {mindist=distance; nearest=predicted[j].id;}
      }
      observations[i].id = nearest;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   default_random_engine gen;
   double w_normalizer=0.0;

   for (unsigned int i=0; i<num_particles;++i)
   {
     // the particle x, y coordinates
     double p_x = particles[i].x;
     double p_y = particles[i].y;
     double p_theta = particles[i].theta;

     // A vector to hold the map landmark locations predicted to be within sensor range of the particle
    vector<LandmarkObs> predictions;

    // for each map landmark
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {

         // get id and x,y coordinates
         double lm_x = map_landmarks.landmark_list[j].x_f;
         double lm_y = map_landmarks.landmark_list[j].y_f;
         int lm_id = map_landmarks.landmark_list[j].id_i;

         //the landmark is in the sensor range.
         double distance = dist(lm_x, lm_y, p_x, p_y);
         if (distance <= sensor_range) {
           predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
         }
      }

       // observation vector
       vector<LandmarkObs> t_obs;
       for (unsigned int j = 0; j < observations.size(); j++) {
         //transformation to map coordinates
         double idd=observations[j].id;
         double map_x = cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y + p_x;
         double map_y = sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y + p_y;
         t_obs.push_back(LandmarkObs{idd, map_x, map_y});
       }

       dataAssociation(predictions, t_obs);

       particles[i].weight = 1.0;

       double o_x, o_y, pr_x, pr_y;

       for (unsigned int j = 0; j < t_obs.size(); ++j)
       {

         for (unsigned int k = 0; k < predictions.size(); k++) {
           if (predictions[k].id == t_obs[j].id) {
             pr_x = predictions[k].x;
             pr_y = predictions[k].y;
           }
         }

         o_x = t_obs[j].x;
         o_y = t_obs[j].y;
         // calculate weight using multivariate Gaussian
         double std_x = std_landmark[0];
         double std_y = std_landmark[1];
         double weight =
(1/(2*pi*std_x*std_y)) * exp(-((pr_x-o_x)*(pr_x-o_x)/(2*(std_x)*(std_x)) + (pr_y-o_y)*(pr_y-o_y)/(2*(std_y)*(std_y))));

         particles[i].weight *= weight;
       }

       //w_normalizer+= particles[i].weight;
     }

     //to normailze updateWeights
     //for (unsigned int i = 0; i < particles.size(); ++i) {
      //  particles[i].weight /= w_normalizer;
    //  }
  }

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
default_random_engine gen;
     vector<Particle> new_particles;
     vector<double> weights;
     for (unsigned int i = 0; i < num_particles; ++i) {
       weights.push_back(particles[i].weight);
     }

   	//random particle index
   	uniform_int_distribution<int> particle_index(0, num_particles - 1);
   	int current_index = particle_index(gen);

   	double beta = 0.0;
    double max_weight = *max_element(weights.begin(), weights.end());
    uniform_real_distribution<double> random_weight(0.0, max_weight);

   	for (int i = 0; i < num_particles; ++i) {
   		beta += random_weight(gen)*2.0;
   	  while (beta > weights[current_index]) {
   	    beta -= weights[current_index];
   	    current_index = (current_index + 1) % num_particles;
   	  }
   	  new_particles.push_back(particles[current_index]);
   	}
   	particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
