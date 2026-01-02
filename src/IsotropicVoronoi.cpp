// CM: Updated by Mithila.
// Copyright (C) 2018, Keith Pedersen (Keith.David.Pedersen@gmail.com)
// All rights reserved

// This file defines some useful worker functions for creating isotropic vectors
// and using MonteCarlo integration to determine Voronoi area. 

#include "kdpVectors.hpp"

#include "pqRand.hpp"
#include "distributions.hpp"

#include <cmath>
#include <future>
#include <algorithm> // std::min_element

// CM: Specific function for Vec3 vector addition
void vector_add_vec3(std::vector<kdp::Vec3>& a, const std::vector<kdp::Vec3>& b) {
    assert(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i) {
        a[i].x1 += b[i].x1;
        a[i].x2 += b[i].x2;
        a[i].x3 += b[i].x3;
    }
}

// CM: Specific function for size_t vector addition
void vector_add_size_t(std::vector<size_t>& a, const std::vector<size_t>& b) {
    assert(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
}

// CM: Specialized vector addition for size_t
template<>
void vector_add<size_t>(std::vector<size_t>& a, const std::vector<size_t>& b) {
    assert(a.size() == b.size());
    for(size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
}

// Generate an isotropic unit vector via rejection sampling
kdp::Vec3 IsoVec3(pqRand::engine& gen)
{
	// The easiest way to draw an isotropic 3-vector is rejection sampling.
	// It provides slightly better precision than a theta-phi implementation,
	// and is not much slower.
	
	using vec3_t = kdp::Vec3;
	
	vec3_t iso(false); // false means don't initialize the Vec3
	double r2;
	
	do
	{
		// U_even has enough precision for this application
		// because we scale by r2
		iso.x1 = gen.U_even();
		iso.x2 = gen.U_even();
		iso.x3 = gen.U_even();
		
		r2 = iso.Mag2();
	}
	while(r2 > 1.); // Draw only from the unit sphere (in the positive octant)
	
	// It we wanted to scale the length from some distribution, 
	// the current length r**2 is a random variate (an extra DOF since 
	// only two DOF are required for isotropy) which can be 
	// easily converted into a random U(0, 1) via its CDF 
	// 	u = CDF(r2) = (r2)**(3/2)
	// (you can work this out from the differential volume dV/(V * d(r**2))) = 3 sqrt(r2) / 2.
	
	// This U(0,1) can then be plugged into any quantile function Q(u)
	// for the new length. One should check to see if 
	// dividing out by the current length can be rolled into Q(u), 
	// so that you that you use u to draw the scale that the takes 
	// iso to its new length.
	
	// In this application, we simply want unit vectors, 
	// so we throw away the random entropy of r2.
	iso /= std::sqrt(r2);
		
	// Use random signs to move to the other 7 octants
	gen.ApplyRandomSign(iso.x1);
	gen.ApplyRandomSign(iso.x2);
	gen.ApplyRandomSign(iso.x3);
	
	return iso;
}

////////////////////////////////////////////////////////////////////////

static constexpr size_t INVALID = size_t(-1);

// Brute force search for nearest vector; used for validation
// Several monotonic metrics can be used; only one is the fastest (least FLOPS)
size_t NearestNeighbor(kdp::Vec3 const& vec, std::vector<kdp::Vec3> const& candidates)
{
	size_t nearestNeighbor = INVALID;
	{
		double bestMetric = -INFINITY; // INFINITY
		
		for(size_t i = 0; i < candidates.size(); ++i)
		{
			//~ double const metric = vec.InteriorAngle(candidates[i]);
			//~ double const metric = (vec - candidates[i]).Mag2(); // Twice as fast as interior angle
			double const metric = vec.Dot(candidates[i]); // Slightly faster than chord length
			
			if(metric > bestMetric) // (metric < bestMetric)
			{
				nearestNeighbor = i;
				bestMetric = metric;
			}
		}
		assert(nearestNeighbor < candidates.size());
	}
	
	return nearestNeighbor;
}

////////////////////////////////////////////////////////////////////////

// Pass a Voronoi area job to a thread. This class is useful if
// some setup work should be accomplished once (instead of by each thread).
// Multiple threads can share a job.
// A minimal class for passing jobs to threads; a useful paradigm to test several algorithms
// When classes extend this class, they can add additional work to Setup_Extra.
// NOTE: In general, we may want to delay setup till after construction. 
// This way, we can asynchronously launch threads and get them doing their spawning work, 
// then do Setup while those threads are still preparing to work.
// NOTE: Perhaps a better method has the first thread ready do the Setup, 
// or somehow distribute it among threads (probably not practical, but who knows).
// No time to try this now.
class ThreadJob
{
	// Threads directly access these elements ONLY for read (WAAH ... we're all adults here)
	// Use API to edit data
	public:
		std::vector<kdp::Vec3> const* sample;
		size_t sampleSize_Voronoi; // Voronoi area sample size (per thread)

		std::vector<kdp::Vec3> center;
		std::vector<size_t> count;
		
	protected:			
		mutable std::mutex threadLock; // Maintain thread safety
		mutable std::condition_variable isValid; // Threads wait till the job is valid
		
		bool valid; // Is the job valid?
		
		virtual void Setup_Extra() {return;} // Leave room for extra setup
		
	public:
		ThreadJob(): // An empty job
			sample(nullptr), sampleSize_Voronoi(0), valid(false) 
		{return;}
		
		ThreadJob(std::vector<kdp::Vec3> const& sample_in, 
			size_t const sampleSize_Voronoi_in):
		ThreadJob()
		{
			Setup(sample_in, sampleSize_Voronoi_in);
		}
		
  // CM: Updated Accumulate to use specific vector addition functions
  void ThreadJob::Accumulate(std::vector<kdp::Vec3> const& center_plus, 
			     std::vector<size_t> const& count_plus)
  {
    std::lock_guard<std::mutex> lock(threadLock);
    vector_add_vec3(center, center_plus);  // Use the specialized Vec3 addition
    vector_add_size_t(count, count_plus);  // Use the specialized size_t addition
  }
  
  // Wait until the job is valid
		void WaitTillValid() const
		{
			std::unique_lock<std::mutex> lock(threadLock);
			
			while(not valid)
				isValid.wait(lock);
		}
		
		// Initialize the job so the threads can use it
		void Setup(std::vector<kdp::Vec3> const& sample_in, 
			size_t const sampleSize_Voronoi_in)
		{
			{
				std::lock_guard<std::mutex> lock(threadLock);
				valid = false;
				
				sample = &sample_in;
				sampleSize_Voronoi = sampleSize_Voronoi_in;
				
				// Each vector is one sample towards its geometric center and Voronoi area
				center = *sample;
				count.assign(center.size(), 1);
						
				this->Setup_Extra();
				
				valid = true;
			}// Release lock before notifying; no hurry-up-and-wait.
			isValid.notify_all();
		}
};

////////////////////////////////////////////////////////////////////////

// Find the Voronoi area of sample via MonteCarlo integration. 
// Do one piece from the job. Do brute force nearest neighbor search.
void Voronoi_CenterCount(ThreadJob& job, std::string const& seed = "")
{
	job.WaitTillValid();
	
	// Bind reference to sample, tell compiler the size_t is actually constant
	std::vector<kdp::Vec3> const& sample = *(job.sample);
	size_t const sampleSize_Voronoi = job.sampleSize_Voronoi;
	
	if(sample.size())
	{
		// Don't include self count; only done once (not by each thread)
		std::vector<kdp::Vec3> center(sample.size(), kdp::Vec3());
		std::vector<size_t> count(sample.size(), 0);	
		
		pqRand::engine gen(false); // Initialize without seeding
		
		if(seed.size())
			gen.Seed_FromString(seed);
		else
			gen.Seed(); // Seed from system entropy
		
		// For each variate iso, find its nearest neighbor.
		// This means iso landed in that vec's Voronoi area. 
		// Keep a running sum of these for each vec, 
		// to define the geometric center of it's Voronoi area, 
		// and a count, which is proportional to its Voronoi area.
		for(size_t j = 0; j < sampleSize_Voronoi; ++j)
		{
			kdp::Vec3 iso = IsoVec3(gen);
			
			size_t const nearestNeighbor = NearestNeighbor(iso, sample);
			
			center[nearestNeighbor] += iso;
			++count[nearestNeighbor];
		}
		
		job.Accumulate(center, count);	
	}
}

////////////////////////////////////////////////////////////////////////

// Do a Voronoi are search by remembering each particle's nearest neighbor. 
// A detailed description can be found in IsotropicVoronoi.pdf
// We now have an additional task to do in Setup_Extra; 
// for each vector, find and cache a list of its nearest neighbors
class ThreadJob_Memory : public ThreadJob
{
	public:
		// For each vector in sample, it's nearest neighbors.
		std::vector<std::vector<size_t>> nearestNeighbors;
		double R;
		
		ThreadJob_Memory(): ThreadJob() {return;}
		
	protected:	
		virtual void Setup_Extra()
		{
			nearestNeighbors.assign(sample->size(), std::vector<size_t>());
			
			// per the AD, this is the optimal search radius
			R = 2.*std::sqrt(2.) * std::pow(double(sample->size()), -0.25);
			double const cos_R = std::cos(R);
				
			for(size_t i = 0; i < nearestNeighbors.size(); ++i)
			{
				// This array of nearest neighbors is symmetric, does not include self
				for(size_t j = 0; j < i; ++j)
				{	
					if((*sample)[i].Dot((*sample)[j]) > cos_R)
					{
						nearestNeighbors[i].push_back(j);
						nearestNeighbors[j].push_back(i);
					}
				}
			}
			
			// Sort the list for easy skipping of neighbors we've already checked
			for(size_t i = 0; i < nearestNeighbors.size(); ++i)
				std::sort(nearestNeighbors[i].begin(), nearestNeighbors[i].end());
		}
};

////////////////////////////////////////////////////////////////////////

// More efficient than brute force search (which should be reviewed for 
// comments regarding the general procedure here).
// See the accompanying document (AD) for math explanations. 
void Voronoi_CenterCount_Memory(ThreadJob_Memory& job, std::string const& seed = "")
{
	job.WaitTillValid();
	
	std::vector<kdp::Vec3> const& sample = *(job.sample);
	std::vector<std::vector<size_t>> const& nearestNeighbors = job.nearestNeighbors;
	size_t const sampleSize_Voronoi = job.sampleSize_Voronoi;	
	
	if(sample.size())
	{
		// Don't include self; only done once (not by each thread)
		std::vector<kdp::Vec3> center(sample.size(), kdp::Vec3());
		std::vector<size_t> count(sample.size(), 0);		
		
		pqRand::engine gen(false); // Don't seed yet
		
		if(seed.size())
			gen.Seed_FromString(seed);
		else
			gen.Seed(); // Seed from system entropy
			
		// per the AD, this is the optimal search radius
		double const cos_HalfR = std::cos(0.5*job.R);
			
		for(size_t j = 0; j < sampleSize_Voronoi; ++j)
		{
			kdp::Vec3 iso = IsoVec3(gen);
			
			size_t closest_i = INVALID;
			double metric;
			double bestMetric = -INFINITY;
			
			for(size_t i = 0; i < sample.size(); ++i)
			{
				metric = iso.Dot(sample[i]);
				
				// Even if we never enter the following if, this will find the NN via brute force (fail-safe)
				if(metric > bestMetric)
				{
					closest_i = i;
					bestMetric = metric;
				}
																	
				// If iso is within R/2 of sample[i], 
				// then all vectors within R/2 of iso are guaranteed to 
				// be contained in sample[i]'s nearestNeighbor list. 
				if(metric > cos_HalfR)
				{
					if(nearestNeighbors[closest_i].size())
					{
						auto itNeighbor = nearestNeighbors[closest_i].begin();
						auto const itEnd = nearestNeighbors[closest_i].end();
						
						// We've already tested index i and below
						if(nearestNeighbors[closest_i].back() > i)
						{
							// Advance past i; can't run out of bounds, last neighbor > i
							while(*itNeighbor <= i) ++itNeighbor;
							
							//~ printf("%lu\n", i + (itEnd - itNeighbor));
						
							for(; itNeighbor not_eq itEnd; ++itNeighbor)
							{
								metric = iso.Dot(sample[*itNeighbor]);
										
								if(metric > bestMetric)
								{
									closest_i = *itNeighbor;
									bestMetric = metric;
								}
							}
						}
						//~ else printf("%lu\n", i);
					}
							
					// We just searched everything within R/2. 
					// At the very least, sample[i] was within that circle, 
					// and it's the nearest neighbor.
					break;
				}
			}
			assert(closest_i < sample.size());
			assert(bestMetric > cos_HalfR);
			//~ assert(closest_i == NearestNeighbor(iso, sample)); // Validate through brute force
			
			// Add the integration vector to isoVec[min_i], to find the geometric center of sample[i]
			center[closest_i] += iso;
			++(count[closest_i]);
		}
		
		job.Accumulate(center, count);
	}
}

////////////////////////////////////////////////////////////////////////

// Used by old map search (look for "Old intersection code"; this was taking all the time
template<typename T>
std::vector<T> ToSortedSet(std::vector<T> vec)
{
	std::sort(vec.begin(), vec.end());
	
	// std::unique only removes duplicate ADJACENT elements; must sort first
	vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
		
	return vec;
}

////////////////////////////////////////////////////////////////////////

// For the map based nearest neighbor search, we have to create a cache of sample z values. 
// A detailed description can be found in IsotropicVoronoi.pdf.
// (sorted by z) containing also the phi and the vector's identity
class ThreadJob_Map : public ThreadJob
{
	public:
		// We map z positions to a vector's identity (index in list) and phi location
		struct zCache
		{
			double z;
			double phi;
			std::vector<kdp::Vec3>::const_iterator vec;
			
			// we sort only by z value (std::sort)
			bool operator< (zCache const& that) const
			{
				return (this->z < that.z);
			}
			
			// std::lower_bound requires bi-directional comparison to double
			bool operator< (double const thatZ) const
			{
				return (this->z < thatZ);
			}
			
			// This is a friend operator (not a member of the class)
			friend bool operator < (double const thisZ, zCache const& that)
			{
				return thisZ < that.z;
			}
		};
		
		std::vector<zCache> zMap; // WAAH: use only for read
		
		ThreadJob_Map(): ThreadJob() {return;}
		
		inline static bool InRange(double const val, std::pair<double, double> const& range)
		{
			// This is better; (not OutOfRange) will say that nan is in range
			return ((val >= range.first) and (val <= range.second));
		}
		
		inline static bool OutOfRange(double const val, std::pair<double, double> const& range)
		{
			return (val < range.first) or (val > range.second);
		}
		
	protected:
	
		virtual void Setup_Extra()
		{
			zMap.reserve(sample->size() + 2);
			zMap.clear();
			
			// Protect iterating out of bounds or special edge handling
			// by padding the z values with too large values
			zMap.push_back(zCache{-2., INFINITY, sample->cend()});
					
			for(auto it = sample->cbegin(); it not_eq sample->cend(); ++it)
				zMap.push_back(zCache{it->x3, it->Phi(), it});
			
			zMap.push_back(zCache{2., INFINITY, sample->cend()});
			
			std::sort(zMap.begin(), zMap.end());
		}
};

////////////////////////////////////////////////////////////////////////

enum class RegionType : size_t {Polar = 0, MeridianWrap = 1, Standard = 2};

// A more efficient method than memory; search a binary map for nearest neighbors.
// Still uses memory, but scales proportional to sample.size()
void Voronoi_CenterCount_Map(ThreadJob_Map& job, std::string const& seed = "")
{
	job.WaitTillValid();
	
	std::vector<kdp::Vec3> const& sample = *(job.sample);
	std::vector<ThreadJob_Map::zCache> const& zMap = job.zMap;
	size_t const sampleSize_Voronoi = job.sampleSize_Voronoi;
	
	if(sample.size())
	{
		// Don't include self; only done once (not by each thread)
		std::vector<kdp::Vec3> center(sample.size(), kdp::Vec3());
		std::vector<size_t> count(sample.size(), 0);		
		
		pqRand::engine gen(false); // Don't seed yet
		
		if(seed.size())
			gen.Seed_FromString(seed);
		else
			gen.Seed(); // Seed from system entropy
		
		auto INVALID_IT = sample.cend();

		for(size_t j = 0; j < sampleSize_Voronoi; ++j)
		{
			kdp::Vec3 iso = IsoVec3(gen);
			auto nearestNeighbor = INVALID_IT;
			
			{
				double const& z = iso.x3; // iso is local, reference is better
				double const phi = iso.Phi();
				double const sinSquaredTheta = iso.T().Mag2();
				
				// A search region can be defined by z and phi boundaries
				std::pair<double, double> boundary_z, boundary_phi;
				
				// The circular search region is defined by its fraction of the unit sphere (A / (4 Pi)).
				// This choice finds nearly every vector in the first try (for isotropic sample)
				// Notice that it is immediately doubled in the first iteration.
				double searchFraction = 1./double(sample.size());
							
				do
				{
					searchFraction *= 2.; // expand the search each iteration
										
					//~ double const R = 2.*std::asin(std::sqrt(searchFraction));
					// sin(2*x) = 2*sin(x)cos(x)
					// cos(2*x) = cos(x)**2 - sin(x)**2						
					double const cosR = 1. - 2.*searchFraction;
					double const sinSquaredR = 4. * searchFraction * (1. - searchFraction);
										
					// Determine the RegionType and set the phi boundaries
					RegionType const regionType = [&]()
					{
						if(sinSquaredTheta > sinSquaredR)
						{
							double const deltaPhi = 2.*std::asin(std::sqrt(sinSquaredR / 
								(2. * (sinSquaredTheta  + std::sqrt(sinSquaredTheta) * 
								std::sqrt(sinSquaredTheta - sinSquaredR)))));
							
							boundary_phi.first = phi - deltaPhi;
							boundary_phi.second = phi + deltaPhi;
							
							// Detecting a search region wrapping around the far meridian
							// In this case, we wrap the offending phi into [-pi, pi],
							// so that the boundaries now define the opposite of the search regions
							if(boundary_phi.first < -M_PI)
							{
								boundary_phi.first += 2.*M_PI;
								return RegionType::MeridianWrap;
							}
							else if(boundary_phi.second > M_PI)
							{
								boundary_phi.second -= 2.*M_PI;
								return RegionType::MeridianWrap;
							}
							else return RegionType::Standard;
						}
						else return RegionType::Polar;
						// When the pole is in the search region, all phi are valid, 
						// so we don't set or use the phi boundaries
					}();
					
					// z boundaries
					{
						double const partA = z * cosR;
						double const partB = std::sqrt(sinSquaredTheta * sinSquaredR);
						
						boundary_z.first = partA - partB;
						boundary_z.second = partA + partB;
						
						assert(boundary_z.first >= -1.);
						assert(boundary_z.second <= 1.);
						
						// When the pole is included, take everything outside the more central z-boundary
						if(regionType == RegionType::Polar)
						{
							if(z > 0.)
								boundary_z.second = 1.;
							else
								boundary_z.first = -1.;
						}
					}
					
					double bestMetric = cosR; // Only accept candidates from the search region
					{
						// Find the first candidate not less than boundary_z.first
						// We assume that this efficiently chooses between linear and binary search.
						auto itCand = std::lower_bound(zMap.begin(), zMap.end(), 
							boundary_z.first);
						
						// Ensure we start somewhere real
						assert(itCand not_eq zMap.end());
						// Ensure that iterating without an explicit bounds check is safe
						assert(boundary_z.second < zMap.back().z);
						
						for(; itCand->z < boundary_z.second; ++itCand)
						{
							// If the pole is not in the search region, we can
							// exclude the candidate by its phi value.
							// To exclude, we "continue" to skip the metric calculation
							// These jumps are highly predictable ... same for all candidates
							switch(regionType)
							{
								case RegionType::Polar: // do nothing
								break;
								
								case RegionType::MeridianWrap:
									if(not ThreadJob_Map::OutOfRange(itCand->phi, boundary_phi))
										continue;
								break;
								
								case RegionType::Standard:
									if(not ThreadJob_Map::InRange(itCand->phi, boundary_phi))
										continue;
								break;
							}
							
							double const metric = itCand->vec->Dot(iso);
						
							if(metric > bestMetric)
							{
								nearestNeighbor = itCand->vec;
								bestMetric = metric;
							}
						}
					}
				}while(nearestNeighbor == INVALID_IT);
			}
			
			size_t const nearestNeighbor_index = nearestNeighbor - sample.cbegin();
			assert(nearestNeighbor_index < job.sample->size());
			
			// Validate with brute force search	
			//~ assert(nearestNeighbor_index == NearestNeighbor(iso, sample));
			
			center[nearestNeighbor_index] += iso;
			++count[nearestNeighbor_index];
		}// End Monte Carlo Voronoi sampling
		
		job.Accumulate(center, count);
	}
}

// Old set intersection code 
// Sort the lists and remove duplicates ... creating a set
//~ candidates_z = ToSortedSet(std::move(candidates_z));

//~ if(candidates_phi.empty()) // Every phi is valid; so take all z candidates
	//~ candidates = std::move(candidates_z); 
//~ else // Valid candidates must appear in both lists
//~ {
	//~ candidates_phi = ToSortedSet(std::move(candidates_phi));
	
	//~ // set_intersection() writes the intersection to a pre-allocated list;
	//~ // make enough room for the maximal intersection (using invalid indices)
	//~ candidates.assign(std::min(candidates_phi.size(), candidates_z.size()), INVALID);
	
	//~ // Form the intersection and resize it to the correct size
	//~ candidates.resize(
		//~ std::set_intersection(candidates_phi.begin(), candidates_phi.end(),
			//~ candidates_z.begin(), candidates_z.end(), 
			//~ candidates.begin()) - candidates.begin());
////////////////////////////////////////////////////////////////////////

// Map a set sample to the CM frame (Note: copy the incoming vector, so we can alter it)
std::vector<kdp::Vec3> MapToCM(std::vector<kdp::Vec3> sample, double const finalWeight = 1.)
{
	kdp::Vec3 totalVec;
	std::vector<double> lengths;
	lengths.reserve(sample.size());
	
	// Add up the vector total, and cache the lengths
	for(auto const& vec : sample)
	{
		totalVec += vec;
		lengths.emplace_back(vec.Mag());
	}
	
	{
		double const totalWeight = kdp::BinaryAccumulate(lengths);
		
		// Subtract total from each vector, proportional to its length; 
		// this maps the sample to its CM frame in an egalitarian way.
		for(size_t i = 0; i < sample.size(); ++i)
		{
			sample[i] -= totalVec * (lengths[i] / totalWeight);
			lengths[i] = sample[i].Mag(); // Store the new length
		}
	}
	
	double const scaleToFinalWeight = finalWeight / kdp::BinaryAccumulate_Destructive(lengths);
	
	// Re-scale the vectors so they have the total weight requested
	for(auto& vec : sample)
		vec *= scaleToFinalWeight;
	
	return sample;
}

/////////////////////////-///////////////////////////////////////////////

// Draw isotropic Voronoi vectors
// 1. Randomly choose isotropic unit vectors (rawIsoVec)
// 2. Use Monte Carlo integration (iso Vec3) to find Voronoi area of each rawIsoVec
// 3. Correct the vectors the CM frame (explained within)
std::vector<kdp::Vec3> MakeIsoVoronoiVectors(pqRand::engine& gen, 
	size_t const sampleSize, size_t const sampleSize_area, bool const fix_CM = true)
{
	ThreadJob* job = nullptr;
	
	{
		std::vector<kdp::Vec3> rawIsoVec;
		rawIsoVec.reserve(sampleSize);
			
		// Draw random, isotropic vectors --- these are the sampled particles
		// We will then find the Voronoi cells which surround each
		while(rawIsoVec.size() < sampleSize)
			rawIsoVec.push_back(IsoVec3(gen));
			
		// Small job; not worth the extra effort
		// (also, the map search is not fully validated for very large R)		
		if(sampleSize < 20) 
		{
			job = new ThreadJob(rawIsoVec, sampleSize_area);
			
			Voronoi_CenterCount(std::ref(*job), gen.GetState());
		}
		else
		{
			using T = ThreadJob_Map;			
			
			job = new T(); // Wait till threads are launched and spinning up before setting up.
			// Perhaps a better method has the first thread ready do the setup work
			T& myJob = *static_cast<T*>(job);
			
			std::vector<std::thread> threads;
			static constexpr size_t numThreads = 4; // MUST be power of 2 to not lose samples
			
			for(size_t i = 0; i < numThreads; ++i)
			{
				threads.emplace_back(&Voronoi_CenterCount_Map, 
					std::ref(myJob), gen.GetState());
				gen.Jump();
			}
							
			myJob.Setup(rawIsoVec, sampleSize_area / numThreads);
			
			for(size_t i = 0; i < numThreads; ++i)
				threads[i].join();
		}
	}
	
	{
		// The length of centerCount.first is made proportional to its areas
		double const totalCount = double(sampleSize_area + sampleSize);
		
		for(size_t i = 0; i < sampleSize; ++i)
			job->center[i] *= double(job->count[i])/ (totalCount *  job->center[i].Mag());
	}		
	
	if(fix_CM)
		return MapToCM(std::move(job->center));
	else
		return job->center;
}

////////////////////////////////////////////////////////////////////////

/* Create an "isotropic" partition of the sphere,
 * with two polar caps of polar width dTheta/2
 * and an integer number of azimuthal belts with polar width dTheta
 * 
 * 0      <-dTheta->                            Pi
 * |_____|__________|__________|__________|_____|
 *                      theta
 * 
 * dTheta = Pi / (nBelta + 1)
 * 
 * Solid angle of arbitrary tile: Omega = 2 dPhi sin(0.5 dPhi) sin(theta_center)
 * Solid angle of polar cap:      Omega0 = 4 Pi sin**2(dTheta / 4)
 * 
 * For each band, solve dPhi so that Omega ~= Omega0, 
 * then adjust dPhi so their are an integer number of tiles. * 
 */ 
std::vector<std::pair<std::vector<kdp::Vec3>, double>> MakeIsoPartition(size_t const nTarget)
{
	// Given some target number of tile in the partition, 
	// we must find an integer number of bands. 
	size_t const nBelts = [&]()
	{
		// So we solve for the dTheta which gives us 4 Pi / Omega0 = nTarget
		double const dTheta = 4.*std::asin(std::sqrt(1./double(nTarget)));
		// Given this dTheta, round to an integer number of belts
		assert(dTheta > M_PI);
		return size_t(std::ceil(M_PI / dTheta - 1.));
	}();
	double const dTheta = M_PI / double(nBelts + 1);
	
	// This shows up repeatedly when we calculate nPhi
	// dPhi = Omega0 / (2 sin(0.5 * dTheta) sin(theta));
	// nCells = 2Pi / dPhi = (sin(0.5 * dTheta) * sin(theta) / sin**2(0.25 * dTheta)
	double const nPhi_factor = std::sin(0.5*dTheta) / 
		kdp::Squared(std::sin(0.25*dTheta));
	
	using Vec3 = kdp::Vec3;
	std::vector<std::pair<std::vector<kdp::Vec3>, double>> tileVec;
	
	// Emplace the polar caps
	tileVec.emplace_back(std::vector<kdp::Vec3>(), 4.*M_PI*kdp::Squared(std::sin(0.25 * dTheta)));
	tileVec.back().first = {Vec3(0., 0., 1.), Vec3(0., 0., -1.)};
		
	for(size_t iBelt = 0; iBelt < nBelts; ++iBelt)
	{
		double const theta = dTheta * double(iBelt + 1); // central theta
		size_t const nPhi = size_t(std::round(nPhi_factor * std::sin(theta)));
		double const dPhi = 2.*(M_PI / double(nPhi));
		double const Omega = 2.*dPhi * std::sin(0.5*dTheta) * std::sin(theta);
		
		// Create an offset to the phi origin
		double const phi0 = -M_PI + (bool(iBelt & 1) ? 0. : 0.5*dPhi);
		
		tileVec.emplace_back(std::vector<Vec3>(), Omega);
		
		for(size_t k = 0; k < nPhi; ++k)
			tileVec.back().first.emplace_back(1., theta, phi0 + double(k)*dPhi, kdp::Vec3from::LengthThetaPhi);		
	}
	
	// Verify solid angles all approximately the same, right number of vectors
	//~ size_t count = 0;
	
	//~ for(auto const& tile : tileVec)
	//~ {
		//~ printf("%.3e\n", tile.second);
		//~ count += tile.first.size();
	//~ }
	//~ printf("%lu\n", count);
	
	return tileVec;
}

////////////////////////////////////////////////////////////////////////
// Voronoi area rejection sampling -- works well
double constexpr b_Voronoi = 3.592; // Empirical constant from thesis

double f(double const x) // Voronoi area PDF
{	
	static double const log_bFactor = b_Voronoi * std::log(b_Voronoi) - std::lgamma(b_Voronoi);
	
	return std::exp( -b_Voronoi*x + (b_Voronoi - 1.)*std::log(x) + log_bFactor);
}

std::vector<double> SampleVoronoiF(pqRand::engine& gen, size_t const sampleSize)
{
	double constexpr x0_Voronoi = 1.273;
	double constexpr d = (1. + (x0_Voronoi - 1.)*b_Voronoi)/x0_Voronoi;

	double const f_x0 = f(x0_Voronoi);
	double const f_mode = f((b_Voronoi - 1.)/b_Voronoi);
	double const tail_integral = f(x0_Voronoi) / (1. + (x0_Voronoi - 1.)*b_Voronoi);
	double const R = f_mode / (f_mode + tail_integral);
	//~ double const R = 0.7416188662047427;
		
	std::vector<double> sample;
	double totalF = 0.;
	
	pqRand::exponential tail(d);
	double x = 1.;
	bool good = false;

	while(sample.size() < sampleSize)
	{
		// If R is the ratio of proposal function areas, 
		// then we must re-choose the region each time we try a variate.
		// If we choose the region, then keep drawing from that region until 
		// we get a good variate, we bias the sample to the proposal distribution.
				
		if(gen.U_even() < R) // mode
		{
			x = gen.U_uneven() * x0_Voronoi;			
			good = (gen.U_uneven() * f_mode < f(x));
		}
		else // tail
		{
			x = tail(gen) + x0_Voronoi;
			good = (gen.U_uneven() * f_x0 * std::exp(-d*(x - x0_Voronoi)) < f(x));
		}
		
		if(good)
		{
			sample.push_back(x);
			totalF += x;
			good = false;
		}
		
		//~ if(gen.U_even() < R) // mode
		//~ {
			//~ do
				//~ x = gen.U_uneven() * x0;
			//~ while(gen.U_uneven() * f_mode > f(x));
			
			//~ sample.push_back(x);
		//~ }
		//~ else // tail
		//~ {
			//~ do
				//~ x = tail(gen);
			//~ while(gen.U_uneven() * f_x0 * std::exp(-d*x) > f(x + x0));
			
			//~ sample.push_back(x + x0);
		//~ }
		
		//~ totalF += sample.back();
	}
	
	// Rescale f
	for(auto& f : sample)
		f /= totalF;
		
	return sample;
}

////////////////////////////////////////////////////////////////////////

// Find each particle's nearest neighbor (brute force), returning a vector of distances and indices.
std::vector<std::pair<double, size_t>> NearestNeighborDistance(std::vector<kdp::Vec3> const& particles)
	//~ std::string const& name = "angles")
{
	// This historgram told us that, in a grid, the same distances appear again and again
	//~ kdp::BinSpecs angle("angle", 100, {0, 0.5});
	//~ kdp::KDPHistogram1 hist(name, angle);
	
	auto distances = std::vector<std::vector<double>>(particles.size(), 
		std::vector<double>(particles.size(), -1.));
	// Initialize to -1 so min_element will catch unfilled 
	
	for(size_t i = 0; i < particles.size(); ++i)
	{
		for(size_t j = 0; j < i; ++j)
		{
			double const r2 = particles[i].InteriorAngle(particles[j]);
			//~ hist.Fill(r2);
		
			// The distance metric is symmetric
			distances[i][j] = r2;
			distances[j][i] = r2;
			
			// However, we read out the smallest value in each ROW, 
			// so we must populate the whole matrix to figure this out.
		}
		distances[i][i] = INFINITY;
	}
	
	std::vector<std::pair<double, size_t>> minDist;
	
	for(size_t i = 0; i < particles.size(); ++i)
	{
		size_t const min_j = std::min_element(distances[i].cbegin(), distances[i].cend()) 
			- distances[i].cbegin();
		minDist.emplace_back(distances[i][min_j], min_j);
	}
	
	//~ std::sort(minDist.begin(), minDist.end());
	
	return minDist;
}
