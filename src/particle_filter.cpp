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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
								// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
								//   x, y, theta and their uncertainties from GPS) and all weights to 1.
								// Add random Gaussian noise to each particle.
								// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
								// std::cout << "in init" << std::endl;
								if (is_initialized) {
																return;
								}
								num_particles = 100;

								double std_x = std[0];
								double std_y = std[1];
								double std_theta = std[2];

								normal_distribution<double> dist_x(x, std_x);
								normal_distribution<double> dist_y(y, std_y);
								normal_distribution<double> dist_theta(theta, std_theta);

								for (int i=0; i<num_particles; i++) {
																Particle particle;

																particle.id = i;
																particle.x = dist_x(gen);
																particle.y = dist_y(gen);
																particle.theta = dist_theta(gen);
																particle.weight = 1.0;

																// particle.x += dist_x(gen);
																// particle.y += dist_y(gen);
																// particle.theta += dist_theta(gen);

																particles.push_back(particle);
								}
								is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
								// TODO: Add measurements to each particle and add random Gaussian noise.
								// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
								//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
								//  http://www.cplusplus.com/reference/random/default_random_engine/
								// std::cout << "in Prediction " << std::endl;
								double std_x = std_pos[0];
								double std_y = std_pos[1];
								double std_theta = std_pos[2];

								normal_distribution<double> dist_x(0, std_x);
								normal_distribution<double> dist_y(0, std_y);
								normal_distribution<double> dist_theta(0, std_theta);

								for (int i=0; i<num_particles; i++) {
																// when yaw_rate is near constant update only x and y values
																if (fabs(yaw_rate) < 0.00001) {
																								particles[i].x += velocity * delta_t * cos(particles[i].theta);
																								particles[i].y += velocity * delta_t * sin(particles[i].theta);
																}
																else {
																								particles[i].x += velocity / yaw_rate * ( sin( particles[i].theta + yaw_rate * delta_t ) - sin( particles[i].theta ) );
																								particles[i].y += velocity / yaw_rate * ( cos( particles[i].theta ) - cos( particles[i].theta + yaw_rate * delta_t ) );
																								particles[i].theta += yaw_rate * delta_t;
																}

																//Add noise
																particles[i].x += dist_x(gen);
																particles[i].y += dist_y(gen);
																particles[i].theta += dist_theta(gen);
								}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
								// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
								//   observed measurement to this particular landmark.
								// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
								//   implement this method and use it as a helper during the updateWeights phase.
								// std::cout << "in DA" << std::endl;
								for (unsigned int i=0; i<observations.size(); i++) {
																double min_distance = numeric_limits<double>::max();
																int map_id = -1;

																for (unsigned int j=0; j<predicted.size(); j++) {
																								double x_distance = observations[i].x - predicted[j].x;
																								double y_distance = observations[i].y - predicted[j].y;
																								double distance = x_distance * x_distance + y_distance * y_distance;

																								if(distance < min_distance) {
																																min_distance = distance;
																																map_id = predicted[j].id;
																								}
																}
																observations[i].id = map_id;
								}
}

std::vector<LandmarkObs> ParticleFilter::inRangeLandmarks(double x, double y, double sensor_range, const Map map_landmarks){
								// std::cout << "in iRL" << std::endl;
								std::vector<LandmarkObs> landmark_objects;
								double total_sensor_range = sensor_range * sensor_range;

								for(unsigned int i=0; i<map_landmarks.landmark_list.size(); i++) {
																float landmark_x = map_landmarks.landmark_list[i].x_f;
																float landmark_y = map_landmarks.landmark_list[i].y_f;
																int id = map_landmarks.landmark_list[i].id_i;
																double dX = x - landmark_x;
																double dY = y - landmark_y;

																if(dX * dX + dY * dY <= total_sensor_range) {
																								landmark_objects.push_back(LandmarkObs{id, landmark_x, landmark_y});
																}
								}
								return landmark_objects;
}

std::vector<LandmarkObs> ParticleFilter::mapObservations(double x, double y, double theta, const std::vector<LandmarkObs> &observations){
								// std::cout << "in MO" << std::endl;
								std::vector<LandmarkObs> mapped_observations;

								for(unsigned int j = 0; j < observations.size(); j++) {
																double landmark_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
																double landmark_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
																mapped_observations.push_back(LandmarkObs{ observations[j].id, landmark_x, landmark_y });
								}
								return mapped_observations;
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
								// std::cout << "in UW" << std::endl;
								double std_landmark_range = std_landmark[0];
								double std_landmark_bearing = std_landmark[1];

								for (int i = 0; i < num_particles; i++) {
																double x = particles[i].x;
																double y = particles[i].y;
																double theta = particles[i].theta;

																std::vector<LandmarkObs> in_range_landmarks = inRangeLandmarks(x, y, sensor_range, map_landmarks);

																std::vector<LandmarkObs> mapped_observations = mapObservations(x, y, theta, observations);

																dataAssociation(in_range_landmarks, mapped_observations);

																particles[i].weight = 1.0;
																for(unsigned int j = 0; j < mapped_observations.size(); j++) {
																								double observation_x = mapped_observations[j].x;
																								double observation_y = mapped_observations[j].y;
																								int landmark_id = mapped_observations[j].id;

																								double landmark_x, landmark_y;
																								unsigned int k = 0;
																								bool found = false;

																								while( !found && k < in_range_landmarks.size() ) {
																																if ( in_range_landmarks[k].id == landmark_id) {
																																								found = true;
																																								landmark_x = in_range_landmarks[k].x;
																																								landmark_y = in_range_landmarks[k].y;
																																}
																																k++;
																								}

																								double dX = observation_x - landmark_x;
																								double dY = observation_y - landmark_y;

																								double weight = ( 1/(2*M_PI*std_landmark_range*std_landmark_bearing)) * exp( -( dX*dX/(2*std_landmark_range*std_landmark_range) + (dY*dY/(2*std_landmark_bearing*std_landmark_bearing)) ) );

																								// std::cout << "weight" << weight << std::endl;
																								// std::cout << "particle weight" << particles[i].weight << std::endl;
																								if (weight == 0) {
																																particles[i].weight *= 0.00001;
																								} else {
																																particles[i].weight *= weight;
																								}
																								// std::cout << "new weight" << particles[i].weight << std::endl;
																}
								}
}

void ParticleFilter::resample() {
								// TODO: Resample particles with replacement with probability proportional to their weight.
								// NOTE: You may find std::discrete_distribution helpful here.
								//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
								// Store all weights in a new vector and get maximum weight
								// std::cout << "in resample" << std::endl;
								std::vector<double> weights;
								double max_weight = numeric_limits<double>::min();
								for(int i = 0; i < num_particles; i++) {
																weights.push_back(particles[i].weight);
																if ( particles[i].weight > max_weight ) {
																								max_weight = particles[i].weight;
																}
								}
								// std::cout << "got max weight" << max_weight << std::endl;
								// for (int i=0; i<weights.size(); i++){
								//  std::cout << weights[i] << " " << std::flush;
								// }

								// Create distributions
								uniform_real_distribution<double> distDouble(0.0, max_weight);
								uniform_int_distribution<int> distInt(0, num_particles - 1);

								// "The wheel" to resample points
								// Pick a random index and set beta to zero
								int index = distInt(gen);
								double beta = 0.0;

								std::vector<Particle> resampledParticles;
								for(int i = 0; i < num_particles; i++) {
																beta += distDouble(gen) * 2.0;
																while( beta > weights[index]) {
																								// std::cout << index << std::endl;
																								// std::cout << beta << "-" << weights[index] << std::endl;
																								beta -= weights[index];
																								index = (index + 1) % num_particles;
																}
																resampledParticles.push_back(particles[index]);
								}

								particles = resampledParticles;
								// std::cout << "returning from resample" << std::endl;
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
								s = s.substr(0, s.length()-1); // get rid of the trailing space
								return s;
}
string ParticleFilter::getSenseX(Particle best)
{
								vector<double> v = best.sense_x;
								stringstream ss;
								copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
								string s = ss.str();
								s = s.substr(0, s.length()-1); // get rid of the trailing space
								return s;
}
string ParticleFilter::getSenseY(Particle best)
{
								vector<double> v = best.sense_y;
								stringstream ss;
								copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
								string s = ss.str();
								s = s.substr(0, s.length()-1); // get rid of the trailing space
								return s;
}
