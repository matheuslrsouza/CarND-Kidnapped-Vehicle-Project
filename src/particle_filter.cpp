/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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


#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		
		particles.push_back(particle);

		weights.push_back(1.0);
	}

	cout << "is inited!!"<<endl;

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;

	for (int i = 0; i < particles.size(); i++) {
		Particle &p = particles[i];

		//avoid division by zero    
    if (fabs(yaw_rate) > 0.0001) {
			p.x = p.x + (velocity / yaw_rate) * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
			p.y = p.y + (velocity / yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
			p.theta = p.theta + yaw_rate * delta_t;
		} else {
			p.x = p.x + velocity * cos(p.theta) * delta_t;
      p.y = p.y + velocity * sin(p.theta) * delta_t;
		}

		normal_distribution<double> dist_x(p.x, std_pos[0]);
		normal_distribution<double> dist_y(p.y, std_pos[1]);
		normal_distribution<double> dist_theta(p.theta, std_pos[2]);

		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int iPart = 0; iPart < particles.size(); iPart++) {

		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;

		Particle &p = particles[iPart];
		double final_weight = 1.0;

		for (int iObs = 0; iObs < observations.size(); iObs++) {
			
			LandmarkObs obs = observations[iObs];
			
			// transform to map x coordinate
			double x_map = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);

			// transform to map y coordinate
			double y_map = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);

			// check which landmark is closest to predicted observation using nearest neighbor
			double best_dist = 1000;
			Map::single_landmark_s best_landmark;
			for (int iLand = 0; iLand < map_landmarks.landmark_list.size(); iLand++) {
				Map::single_landmark_s landmark = map_landmarks.landmark_list[iLand];

				double delta_x = x_map - landmark.x_f;
				double delta_y = y_map - landmark.y_f;

				double dist = sqrt(delta_x*delta_x + delta_y*delta_y);

				if (dist < sensor_range) {
					
					if (dist < best_dist) {
						best_dist = dist;
						best_landmark = landmark;
					}

				}
			}

			// is there association?
			if (best_dist != 1000) {
				associations.push_back(best_landmark.id_i);
				sense_x.push_back(x_map);
				sense_y.push_back(y_map);
			}

			// applying Multivariate_normal_distribution
			// calculate normalization term
			double gauss_norm = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));

			// calculate exponent
			double mu_x = best_landmark.x_f;
			double mu_y = best_landmark.y_f;
			double exponent = (pow(x_map - mu_x, 2))/(2 * pow(std_landmark[0], 2)) + 
											  (pow(y_map - mu_y, 2))/(2 * pow(std_landmark[1], 2));

			// calculate weight using normalization terms and exponent
			final_weight *= gauss_norm * exp(-exponent);
		}

		p.weight = final_weight;
		weights[iPart] = final_weight;

		//cout <<" weight particle ("<<iPart<<"): "<< p.weight <<endl;

		p.associations = associations;
		p.sense_x = sense_x;
		p.sense_y = sense_y;
		
	}

}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// normalizing weights
	vector<double> norm_weights;
	double sum_weights = 0.0;
	for (int i = 0; i < weights.size(); i++) {
		sum_weights += weights[i];
	}

	for (int i = 0; i < weights.size(); i++) {
		norm_weights.push_back(weights[i]/sum_weights);
	}

	weights = norm_weights;

	random_device rd;
  mt19937 gen(rd());
	discrete_distribution<> d(weights.begin(), weights.end());

	std::vector<Particle> new_particles;

	// selects the surviving particles
	for (int i = 0; i < num_particles; i++) {
		int selected_particle = d(gen);
		new_particles.push_back(particles[selected_particle]);
	}
	
	// do the resample
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
