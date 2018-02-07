# Project Introduction
Your robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

In this project you will implement a 2 dimensional particle filter in C++. Your particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step your filter will also get observation and control data. 

## Inputs to the Particle Filter
You can find the inputs to the particle filter in the `data` directory. 

#### The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

## Algorithm

### Init

Creates and initializes the Particles based on estimates of x, y and theta and their uncertainties from GPS. Adds random Gaussian noise to x, y, and theta to each particle. See `ParticleFilter :: init`

### Prediction

Predicts the next position for each particle, adding random Gaussian noise to x, y and theta. See `ParticleFilter::prediction`

### Update

Updates the particle weights based on their relevance to the vehicle position. See `ParticleFilter::updateWeights`

**Transform coordinates**

The observations of the vehicle to the landmarks are sent in the car's coordinate system, where the x is always in the direction to which the vehicle points and y is perpendicular and to the left of x. Knowing this the first thing to do is the transformation of these observations into the map's coordinate system, using the position and direction of the particle. After that step we should execute the landmark association

**Association**

Transformation allows the observation of the vehicle to be applied to each particle. Now that the observations were fitted to the particle, we associate each observation with a landmark (using a nearest-neighbors data association), respecting the range of the sensor, in order not to observe  landmarks out of reach.

**Particle weights**

With each observation associated with a landmark, it is time to measure what weight each particle should receive. We first calculate the probability of each observation using mult-variate Gaussian distribution, where the mean is the position of the associated landmark and x is the particle observation value, the smaller the distance between the observation and the landmark, the greater the likelihood. Then multiply all the calculated measurement probabilities together to obtain the final weight of the particle.

img

### Resample

Using the defined weights, a new collection of particles is created, where the ones with larger weights are more likely to survive than those with smaller weights. The larger particles will probably appear more than once in this new collection.
To facilitate this implementation, the function `std::discrete_distribution` was used. See `ParticleFilter::resample`

