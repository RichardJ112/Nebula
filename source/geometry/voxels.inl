#include "voxels.h"
#include <iostream>
#include <fstream>

namespace nbl { namespace geometry {

template<bool gpu_flag>
inline voxels<gpu_flag>::voxels(real voxel_size, vec3 shape, std::vector<int> initial_geometry, int max_save_height, real sim_depth)
{
	_voxel_size = voxel_size;
	_AABB_min = vec3{ 0, 0, 0 };
	_AABB_max = vec3{ shape.x * voxel_size, shape.y * voxel_size, shape.z * voxel_size };
	
	_size_x = (int)shape.x;
	_size_y = (int)shape.y;
	_size_z = (int)shape.z;

	_max_save_height = max_save_height;
	_min_save_height = max_save_height - 5; // save at least 5 layers of voxels
	
	const vec3 m = _AABB_max - _AABB_min;
	_max_extent = magnitude(m);

	_mat_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_tag_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	_e_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);
	//_dz_grid.resize(int(shape.x) * int(shape.y) * int(shape.z), 0);
	_species_grid.resize(static_cast<int>(shape.x) * static_cast<int>(shape.y) * static_cast<int>(shape.z), 0);

		if (initial_geometry.size() != _size_x * _size_y * _size_z)
		{
			//throw std::invalid_argument("initial geometry of wrong shape");
			std::clog << "Warning: initial geometry of wrong shape\n";
		}
	
	_mat_grid = initial_geometry;

	
}

template<bool gpu_flag>
CPU voxels<gpu_flag> voxels<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	const real VOXEL_SIZE = 0.3; // voxel size in nanometers (0.27 nm is appr. one atom of Si)

	const int SIZE_X = 201; // horizontal size in the x direction in voxels
	const int SIZE_Y = 201; // horizontal size in the y direction in voxels
	const int SIZE_Z = 201; // vertical size in voxels
	const real SIM_DEPTH = 500; // simulation depth under the voxels at z < 0 for SEM bulk samples in nm

	const int SAMPLE_HEIGHT = 101; // height of the sample (length between the sample and the top of the simulation domain, in vacuum) in voxels
	
	vec3 shape = { SIZE_X, SIZE_Y, SIZE_Z };
	
	// Hier moet ergens een grid gemaakt worden, 
		// een voor het materiaal, 
		// een voor the timestamp van de electronen, 
		// een voor de energy van het deposition electron,
		// en een voor de dz van dat electron
	
	 // set the shape 

	// set the initial geometry
	std::vector<int> ini_geom;
	
	ini_geom.resize(SIZE_X * SIZE_Y * SIZE_Z, 0); 
	
	for (int i = 0; i < SIZE_X; i++) {
		for (int j = 0; j < SIZE_Y; j++) {
			for (int k = 0; k < SAMPLE_HEIGHT; k++) {
				ini_geom.at( i + j * SIZE_X + k * SIZE_X * SIZE_Y) = -123;
			}
		}
	}

	// Add layers of detectors for testing
	
	for (int i = 0; i < SIZE_X; i++) {
		for (int j = 0; j < SIZE_Y; j++) {
			ini_geom.at(i + j * SIZE_X + 0 * SIZE_X * SIZE_Y) = -126;
		}
	}
	
	/*for (int i = 0; i < SIZE_X; i++) {
		for (int j = 0; j < SIZE_Y; j++) {
			ini_geom.at(i + j * SIZE_X + (0) * SIZE_X * SIZE_Y) = -126;
		}
	}*/
	
	

	// TODO: error message
	/*if (triangles.empty())
		throw std::runtime_error("No triangles provided!");

	vec3 AABB_min = triangles[0].AABB_min();
	vec3 AABB_max = triangles[0].AABB_max();

	for (const triangle t : triangles)
	{
		const vec3 tri_min = t.AABB_min();
		const vec3 tri_max = t.AABB_max();

		AABB_min =
		{
			std::min(AABB_min.x, tri_min.x),
			std::min(AABB_min.y, tri_min.y),
			std::min(AABB_min.z, tri_min.z)
		};
		AABB_max =
		{
			std::max(AABB_max.x, tri_max.x),
			std::max(AABB_max.y, tri_max.y),
			std::max(AABB_max.z, tri_max.z)
		};
	}

	AABB_min -= vec3{ 1, 1, 1 };
	AABB_max += vec3{ 1, 1, 1 };
	*/
	
	voxels<false> geometry(VOXEL_SIZE, shape, ini_geom, SAMPLE_HEIGHT + 1, SIM_DEPTH);
	return geometry;
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::destroy(voxels<gpu_flag> & geometry)
{
	detail::voxels_factory<gpu_flag>::destroy(geometry);
}

template<bool gpu_flag>
PHYSICS bool voxels<gpu_flag>::in_domain(vec3 pos)
{
	return ((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
		&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
		&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z));
}

template<bool gpu_flag>
PHYSICS intersect_event voxels<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const* ignore_triangle, int ignore_material) const
{
	intersect_event evt{ distance, nullptr };

	// 
	const real x = start.x / _voxel_size; // create vars for the location elements
	const real y = start.y / _voxel_size;
	const real z = start.z / _voxel_size;
	
	const real dx = direction.x; // create vars for the direction elements
	const real dy = direction.y;
	const real dz = direction.z;
	
	//std::clog << "\r                      " << x << "  " << y << "  " << z << "  " << dx << "  " << dy << "  " << dz << "                       ";

	/*if(z + dz*distance > _size_z) // if the electron is under the voxels and it cannot leave it, return the distance-null event
	{
		return evt;
	}*/

	const vec3 dr = { dx, dy, dz };
	vec3 delta_S = { 0, 0, 0 };

	real delta_s_min = distance / _voxel_size; // holds the shortest pathlength to an intersection

	const int start_mat = ignore_material;

	// calculate the distances to the first planes in the 3 dimensions
	real delta_x;	// delta_x is the distance (perpendicular to the plane) to the next plane in the x direction
	if (dx > 0) {
		delta_x = std::ceil(x) - x;
		delta_S.x = delta_x / dx; // delta_S.x now is the path length to the next plane in the x direction
	}
	else if(dx < 0){
		delta_x = x - std::floor(x);
		delta_S.x = -delta_x / dx;
	}
	else // dx == 0
	{
		delta_S.x = distance / _voxel_size;
	}

	real delta_y;
	if (dy > 0) {
		delta_y = std::ceil(y) - y;
		delta_S.y = delta_y / dy;
	}
	else if(dy < 0){
		delta_y = y - std::floor(y);
		delta_S.y = -delta_y / dy;
	}
	else // dy == 0
	{
		delta_S.y = distance / _voxel_size;
	}

	real delta_z;
	if (dz > 0) {
		delta_z = std::ceil(z) - z;
		delta_S.z = delta_z / dz;
	}
	else if(dz < 0) {
		delta_z = z - std::floor(z);
		delta_S.z = -delta_z / dz;
	}
	else // dz == 0
	{
		delta_S.z = distance / _voxel_size;
	}

	if (delta_S.x < 0.000001)
	{
		delta_S.x += std::abs(1 / dx); // add one path length between two planes to delta_S.x
	}
	if (delta_S.y < 0.000001)
	{
		delta_S.y += std::abs(1 / dy);
	}
	if (delta_S.z < 0.000001)
	{
		delta_S.z += std::abs(1 / dz);
	}
	
	
	/*if(z >= _size_z)
	{
		delta_S.z = (z - _size_z) / -dz; // the distance to the z=0 plane is z / (1/dz)

		delta_S.x = std::ceil((((z - _size_z) / -dz) - delta_S.x) * std::abs(dx)) / std::abs(dx) + delta_S.x; // calculate distances to the first intersections with x and y planes 
		delta_S.y = std::ceil((((z - _size_z) / -dz) - delta_S.y) * std::abs(dx)) / std::abs(dy) + delta_S.y; //   such that z > 0

		std::clog << delta_S.x << "   " << delta_S.y << "   " << delta_S.z << "\n";
	}*/


	while (distance / _voxel_size >= delta_s_min) {
		
		//std::clog << "\n" << delta_s_min ;

		// Determine minimum path length from delta_S
		delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z);

		/*if(delta_s_min < 0.000001) // handle the case that delta_s_min is zero, which is typically the case after an intersection
		{
			//std::clog << "\n0 in delta_s_min " << delta_S.x << "  " << delta_S.y << "  " << delta_S.z;
			if(delta_S.x < 0.000001)
			{
				delta_S.x += std::abs(1 / dx); // add one path length between two planes to delta_S.x
			}
			if (delta_S.y < 0.000001)
			{
				delta_S.y += std::abs(1 / dy);
			}
			if (delta_S.z < 0.000001)
			{
				delta_S.z += std::abs(1 / dz);
			}
			delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z); // and determine the minimum again
		}*/

		int min_index = 0; // min_index is 0 for een intersection with the x-plane, 1 for an intersection with the y-plane
		//and 2 for an intersection with the z-plane
		if(delta_s_min == delta_S.y)
		{
			min_index = 1;
		}
		else if (delta_s_min == delta_S.z)
		{
			min_index = 2;
		}
		
		const int min_i = min_index; // store min_index in a constant

		vec3 new_pos = start / _voxel_size + (delta_s_min /*+ 0.001*/) * dr; // new position if there is an intersection in voxels

		vec3 pos = new_pos * _voxel_size; // new position in nm, for check

		if(!((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
			&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
			&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z))) 
			// check whether the electron is still in the simulation domain, if not, return
		{
			//std::clog << "\n exit";
			return evt;
		}

		//std::clog << "\nposition: " << new_pos.x << "  " << new_pos.y << "  " << new_pos.z;

		int k; // indices of the material grid
		int l;
		int m;

		real dx_sgn; // signs of dx, dy and dz
		real dy_sgn;
		real dz_sgn;
		
		//std::clog << "\n" << dx_sgn << "   " << dy_sgn << "   " << dz_sgn;

		// determine the indices of the next cell that the electron will enter
		switch (min_i)
		{
		case 0: // intersection with x-plane
			dx_sgn = dx / std::abs(dx); // determine sign of dx
			k = static_cast<int>(new_pos.x + 0.1 * dx_sgn); // calculate the indices by flooring the new position plus a
			// little inward 
			l = static_cast<int>(new_pos.y);
			m = static_cast<int>(new_pos.z);
			
			if (k >= _size_x || k < 0)
			{
				return evt;
			}
			break;

		case 1: // intersection with y-plane
			dy_sgn =  dy / std::abs(dy);
			k = static_cast<int>(new_pos.x);
			l = static_cast<int>(new_pos.y + 0.1 * dy_sgn);
			m = static_cast<int>(new_pos.z);
			
			if (l >= _size_y || l < 0)
			{
				return evt;
			}
			break;

		default: // intersection with z-plane
			dz_sgn = dz / std::abs(dz);
			k = static_cast<int>(new_pos.x);
			l = static_cast<int>(new_pos.y);
			m = static_cast<int>(new_pos.z + 0.1 *dz_sgn);

			if(m >= _size_z || m < 0)
			{
				return evt;
			}
			break;
		}

		/*if (z > _size_z)
		{
			std::clog << "   (x,y,z): " << x << "   " << y << "   " << z << "\n";
			std::clog << "   (dx,dy,dz): " << dx << "   " << dy << "   " << dz << "\n";
			std::clog << "   (face, dist): " << min_index << "   " << distance / _voxel_size << "\n";
			std::clog << "   new_pos(x,y,z): " << new_pos.x << "   " << new_pos.y << "   " << new_pos.z << "\n";
			std::clog << "   (k,l,m): " << k << "   " << l << "   " << m << "\n";
		}*/
		
		int new_mat = _mat_grid[k + l * _size_x + m * _size_x * _size_y]; // determine material using the material indices

		//std::clog << "   " << new_mat << "   " << start_mat;
		
		if (new_mat != start_mat) { // if thcere is een intersection, return the intersection event

			//std::clog << "\r" << "intersection from " << start_mat << " to " << new_mat << " at " << k << " " << l << " " << m << "    " << _mat_grid[k + l * _size_x + (m-1) * _size_x * _size_y] << " " << _mat_grid[k + l * _size_x + m * _size_x * _size_y] << " " << _mat_grid[k + l * _size_x + (m + 1) * _size_x * _size_y];
			//std::clog << "material at 100 100 299: " << _mat_grid.at(100 + 100 * _size_x + 299 * _size_x * _size_y) << "\n";

			
			evt.isect_distance = (delta_s_min /*+ 0.005*/) * _voxel_size; // set the distance to the intersection

			// Determine voxel side
			int	voxel_side;
			switch (min_i)
			{
			case 0:
				if (dx > 0) {
					voxel_side = 1;
				}
				else {
					voxel_side = 2;
				}
				break;

			case 1:
				if (dy > 0) {
					voxel_side = 3;
				}
				else {
					voxel_side = 4;
				}
				break;

			default:
				if (dz > 0) {
					voxel_side = 5;
				}
				else {
					voxel_side = 6;
				}
				break;
			}

			// put an identifier for the new material and a voxel side identifier in the triangle pointer
			uint64_t isect_id;
			reinterpret_cast<int32_t*>(&isect_id)[0] = new_mat;
			reinterpret_cast<int32_t*>(&isect_id)[1] = voxel_side;
			evt.isect_triangle = reinterpret_cast<triangle*>(isect_id);

			return evt; // return the event
		}

		// if there is no intersection:
		// Calculate new value of delta_S.i for next iteration
		switch (min_i)
		{
		case 0:
			delta_S.x += std::abs(1 / dx);
			break;

		case 1:
			delta_S.y += std::abs(1 / dy);
			break;

		default:
			delta_S.z += std::abs(1 / dz);
			break;
		}
	}
	//std::clog << "\n turn";
	// if no intersection was found, return the standard event. 
	return evt;
	
	/*for (triangle_index_t i = 0; i < _N; ++i)
	{
		if (_triangles + i == ignore_triangle)
			continue;

		// retrieve triangle from global memory
		const triangle this_triangle = _triangles[i];

		int mat_idx_out;
		if (dot_product(this_triangle.get_normal(), direction) < 0)
			mat_idx_out = this_triangle.material_in;
		else
			mat_idx_out = this_triangle.material_out;

		// if the outgoing material is the same as current, nothing happens
		// if the triangle represents a detector which can't see the current
		// particle, nothing happens
		if (mat_idx_out == ignore_material)
			continue;

		const real t = this_triangle.intersect_ray(start, direction);
		if (t > 0 && t < evt.isect_distance)
		{
			evt.isect_distance = t;
			evt.isect_triangle = _triangles + i;
		}
	}*/
	
}

template <bool gpu_flag>
PHYSICS void voxels<gpu_flag>::set_material(vec3 position, int material, int PE_tag, real energy, uint8_t species)
{
	int k = (int) (position.x / _voxel_size);
	int l = (int) (position.y / _voxel_size);
	int m = (int) (position.z / _voxel_size);

	//std::clog << "A deposition at: " << k << "  " << l << "  " << m << "\n";

	//std::clog << _mat_grid[k + l * _size_x + m * _size_x * _size_y] << "   ";

	if(_mat_grid.at(k + l * _size_x + m * _size_x * _size_y) != -123)
	{
		//std::clog << "\n deposition inside material:  " << position.z / _voxel_size;
		return;
	}
	
	_mat_grid.at(k + l * _size_x + m * _size_x * _size_y) = material;
	_tag_grid.at(k + l * _size_x + m * _size_x * _size_y) = PE_tag + 1; // PE_tag begins at 0, so we add 1 to distinguish the deposition from the first electron from the non-deposits
	_e_grid.at(k + l * _size_x + m * _size_x * _size_y) = energy;
	_species_grid.at(k + l * _size_x + m * _size_x * _size_y) = species + 1;

	//std::clog << _mat_grid[k + l * _size_x + m * _size_x * _size_y] << "\n";

	if(m > _max_save_height)
	{
		_max_save_height = m;
	}
	else if (m < _min_save_height)
	{
		_min_save_height = m;
	}
}

template <bool gpu_flag>
int voxels<gpu_flag>::get_material(int position) const
{
	return _mat_grid[position];
}

template <bool gpu_flag>
void voxels<gpu_flag>::deposit(vec3 position, vec3 normal, int material, int PE_tag, real energy, uint8_t species)
{
	vec3 pos = position / _voxel_size + 0.1 * normal;
	
	int k = static_cast<int>(std::floor(pos.x));
	int l = static_cast<int>(std::floor(pos.y));
	int m = static_cast<int>(std::floor(pos.z));

	if (_mat_grid.at(k + l * _size_x + m * _size_x * _size_y) != -123)
	{
		pos = position / _voxel_size - 0.1 * normal;

		k = static_cast<int>(std::floor(pos.x));
		l = static_cast<int>(std::floor(pos.y));
		m = static_cast<int>(std::floor(pos.z));
	}
	if (_mat_grid.at(k + l * _size_x + m * _size_x * _size_y) != -123)
	{
		//std::clog << "\n deposition inside material:  " << _tag_grid.at(k + l * _size_x + m * _size_x * _size_y) << "   " << PE_tag;
		//std::clog << "\n";
		return;
	}

	_mat_grid.at(k + l * _size_x + m * _size_x * _size_y) = material;
	_tag_grid.at(k + l * _size_x + m * _size_x * _size_y) = PE_tag + 1; // PE_tag begins at 0, so we add 1 to distinguish the deposition from the first electron from the non-deposits
	_e_grid.at(k + l * _size_x + m * _size_x * _size_y) = energy;
	_species_grid.at(k + l * _size_x + m * _size_x * _size_y) = species + 1;

	//std::clog << _mat_grid[k + l * _size_x + m * _size_x * _size_y] << "\n";

	if (m > _max_save_height)
	{
		_max_save_height = m;
	}
	else if (m < _min_save_height)
	{
		_min_save_height = m;
	}
}

template<bool gpu_flag>
PHYSICS real voxels<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_min() const
{
	return _AABB_min;
}
template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_max() const
{
	return _AABB_max;
}

template <bool gpu_flag>
void voxels<gpu_flag>::save(const std::string file_name)
{
	std::ofstream file;
	file.open(file_name);
	file << _voxel_size << "\t" << _size_x << "\t" << _size_y << "\t" << _max_save_height - _min_save_height + 1 << "\n";

	for(int i = _min_save_height * _size_x * _size_y; i < (_max_save_height + 1) * _size_x * _size_y; i++)
	{
		file << _mat_grid.at(i) << "\t" << _tag_grid.at(i) << "\t" << _e_grid.at(i) << "\t" << _species_grid.at(i) << "\n";
	}
	
	file.close();
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::set_AABB(vec3 min, vec3 max) 
{
	_AABB_min = min;
	_AABB_max = max;

	const vec3 m = max - min;
	_max_extent = magnitude(m);
}


/*namespace detail
{
	template<>
	struct voxels_factory<false>
	{
		inline static CPU voxels<false> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using voxels_t = voxels<false>;
			using triangle_index_t = voxels_t::triangle_index_t;

			std::vector<int> a;
			voxels_t geometry(3, {5, 5, 6}, a, 9, 9);

			
			/*if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			geometry._triangles = reinterpret_cast<triangle*>(malloc(geometry._N * sizeof(triangle)));
			for (triangle_index_t i = 0; i < triangles.size(); ++i)
			{
				geometry._triangles[i] = triangles[i];
			}

			//geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(voxels<false> & geometry)
		{
			//::free(geometry._triangles);

			//geometry._triangles = nullptr;
			//geometry._N = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<>
	struct voxels_factory<true>
	{
		inline static CPU voxels<true> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using voxels_t = voxels<true>;
			using triangle_index_t = voxels_t::triangle_index_t;

			voxels_t geometry;

			if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			// Copy triangle data to device
			cuda::cuda_new<triangle>(&geometry._triangles, geometry._N);
			cuda::cuda_mem_scope<triangle>(geometry._triangles, geometry._N,
				[&triangles](triangle* device)
			{
				for (triangle_index_t i = 0; i < triangles.size(); ++i)
					device[i] = triangles[i];
			});

			geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(voxels<true> & geometry)
		{
			cudaFree(geometry._triangles);

			geometry._triangles = nullptr;
			geometry._N = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

*/

}} // namespace nbl::geometry


