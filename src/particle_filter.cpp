/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang

 // completed by Reza Kakoee - Jun 24, 2018

 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    	
	default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// setting to 50 particles. it can be tuned
	num_particles = 50;	

	for (int i = 0; i < num_particles; ++i) {
		Particle p;
        p.id = i;//unique id is the index
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
        p.weight = 1;
        particles.push_back(p);
        weights.push_back(1);
	}
    is_initialized = true;

}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

default_random_engine gen;
normal_distribution<double> dist_x(0, std_pos[0]);
normal_distribution<double> dist_y(0, std_pos[1]);
normal_distribution<double> dist_theta(0, std_pos[2]);
for (int i = 0; i < num_particles; ++i) {
		Particle* p;
		p = &particles[i];
		if(abs(yaw_rate)>0.0001){//avoiding div by 0
				
				p->x = p->x + (velocity/yaw_rate)*(sin(p->theta+yaw_rate*delta_t) - sin(p->theta));
				p->y = p->y + (velocity/yaw_rate)*(cos(p->theta)-cos(p->theta+yaw_rate*delta_t));
				p->theta = p->theta + yaw_rate*delta_t;

		}else{

		      p->x += velocity * delta_t * cos(p->theta);
		      p->y += velocity * delta_t * sin(p->theta);
		}
		p->x+=dist_x(gen);
		p->y+=dist_y(gen);
		p->theta+=dist_theta(gen);

	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (unsigned int i= 0;i<observations.size();i++)
	{
		LandmarkObs *obs = &observations[i];
		int n_neighbor=0;
		// for each obs finding the nearest neighbor using Euclidian distance
		double mindist=dist(obs->x,obs->y,predicted[0].x,predicted[0].y);
		for(unsigned int j=1;j<predicted.size();j++){
			double curdist = dist(obs->x,obs->y,predicted[j].x,predicted[j].y);
			if(curdist<mindist){
				mindist=curdist;
				n_neighbor=j;
			}
		}
		//store the index of nearest neigbor for each observation in the obs.id to be used later
		obs->id = n_neighbor;
	}


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// these variables are fixed for all particles. so, calculating once
	double gauss_norm;
	gauss_norm= (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
	double x_denom = 2* (std_landmark[0]*std_landmark[0]);
	double y_denom = 2* (std_landmark[1]*std_landmark[1]);	


	// loop through particles
	for (unsigned int i=0;i<particles.size();i++)
	{
		Particle p = particles[i];
		std::vector<LandmarkObs> landmarks_sensor_range;

 		// building the landmark vector for those landmarks that are in the sensor range
		for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
	        
	        // reduce amount of landmarks to only those in sensor range of the particle
	        double landmark_dist = dist(p.x,p.y, map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
	        if (landmark_dist <= sensor_range) 
	        {
	   	      LandmarkObs lmark;
	   	      lmark.id=k;
	   	      lmark.x = map_landmarks.landmark_list[k].x_f;
	   	      lmark.y = map_landmarks.landmark_list[k].y_f;
			  landmarks_sensor_range.push_back(lmark);
			}
		}
	    // building the observation vector based on the x,y in map coordinates
		std::vector<LandmarkObs> observations_trans;
	    for (unsigned int j = 0; j < observations.size(); ++j) {
	      
	      LandmarkObs lobs;
	      double obs_x, obs_y;
	      // Transforming the observation points 
	      obs_x = p.x+ observations[j].x * cos(p.theta) - observations[j].y * sin(p.theta) ;
	      obs_y = p.y+ observations[j].x * sin(p.theta) + observations[j].y * cos(p.theta) ;
	      lobs.id=j;
	      lobs.x = obs_x;
	      lobs.y = obs_y;
	      observations_trans.push_back(lobs);
	    }
	    // finding nearest landmark for each obs point using the dataAssociation function
   		dataAssociation(landmarks_sensor_range, observations_trans);

   		// calculaitng the weights in mvg
   		double mvg=1;
   		for (unsigned int j = 0; j < observations_trans.size(); ++j) {
   				// the obs.id has the index of nearest neighbor after calling dataAssociation function
		      float m_x = landmarks_sensor_range[observations_trans[j].id].x;
      		  float m_y = landmarks_sensor_range[observations_trans[j].id].y;
      
		      double x_diff = observations_trans[j].x - m_x;
      		  double y_diff = observations_trans[j].y - m_y;
      		  double term = ((x_diff * x_diff) / x_denom) + ((y_diff * y_diff) / y_denom);

      		  mvg = mvg * (gauss_norm * exp(-term)) ;// product of all probabilities 
    	}

    	particles[i].weight = mvg;
    	weights[i] = mvg;
	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Vector for new particles
  vector<Particle> new_particles(num_particles);
  
  // Use discrete distribution to return particles by weight
  default_random_engine gen;
  for (int i = 0; i < num_particles; ++i) {
    discrete_distribution<int> disc(weights.begin(), weights.end());
    new_particles[i] = particles[disc(gen)];
  }
 
  particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
