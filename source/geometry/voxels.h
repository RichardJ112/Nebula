#ifndef __GEOMETRY_VOXELS_H_
#define __GEOMETRY_VOXELS_H_

#include "../core/triangle.h"
#include "../core/events.h"



namespace nbl { namespace geometry {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 */
	template<bool gpu_flag>
	struct voxels_factory;
}

/**
 * \brief Stores geometry as a list of voxel cells.
 *
 * This class is responsible for the collision detection system. It holds the
 * simulation domain (a finite axis-aligned box) and a vector which holds the materials inside the voxels.
 *
 * It is allocated by the static {@link create} and {@link destroy} functions.
 * There is a simple straightforward constructor, but no destructor.
 */
template<bool gpu_flag>
class voxels
{
public:
	using triangle_index_t = uint32_t; ///< Type for indexing the triangles

	/**
	 * \brief Constructor of the voxels class
	 *
	 * \param vox_size Size of a side of the voxels in nanometer
	 * \param shape Shape of the simulation domain in voxels, as vec3 {size_x, size_y, size_z}
	 * \param initial_geometry The initial geometry of the sample, as a std::vector<int> of length size_x*size_y*size_z
	 */
	voxels(real voxel_size, vec3 shape, std::vector<int> initial_geometry, int max_save_height, real sim_depth);
	/**
	 * \brief Allocate memory for the triangles on the correct device.
	 *
	 * \param triangles List of triangles to be used in the simulation.
	 */
	static CPU voxels create(std::vector<triangle> const & triangles);
	/**
	 * \brief Destroy the triangle list, deallocating the data.
	 */
	static CPU void destroy(voxels & geometry);

	/**
	 * \brief Check whether a certain position is part of the simulation domain.
	 *
	 * \param position Position to check.
	 */
	inline PHYSICS bool in_domain(vec3 position);

	/**
	 * \brief Try to move a particle, checking for collisions with triangles.
	 *
	 * Try to propagate from \p start to \p start + \p distance * \p direction.
	 * If there is no collision with a triangle, the returned intersect_event
	 * contains a `nullptr` triangle pointer.
	 *
	 * \param start           The particle's starting position
	 * \param direction       Direction the particle is going in
	 * \param distance        Distance the particle travels, in units of direction
	 * \param ignore_triangle Triangle to ignore. This is not used in this voxel version.
	 * \param ignore_material Destination material to ignore (most likely the
	 *                        particle's current material).
	 *
	 * TODO ignore_material datatype
	 */
	inline PHYSICS intersect_event propagate(vec3 start, vec3 direction, real distance,
		triangle const * ignore_triangle, int ignore_material) const;

	/**
	 * \brief Set the material of a certain voxel
	 *
	 * Set the material in the voxel that contains position to value material. 
	 *
	 * \param position           Where to set the material
	 * \param material			 The material to which the voxel will be set
	 */

	inline PHYSICS void set_material(vec3 position, int material, int PE_tag, real energy, uint8_t species);

	/**
	 * \brief Get the material of a voxel
	 */

	inline PHYSICS int get_material(int position) const;

	/**
	 * \brief Deposit at position position
	 */

	inline PHYSICS void deposit(vec3 position, vec3 normal, int material, int PE_tag, real energy, uint8_t species);

	/**
	 * \brief Get the maximum distance that can be travelled inside the
	 *        simulation domain.
	 */
	inline PHYSICS real get_max_extent() const;

	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_min() const;
	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline PHYSICS vec3 AABB_max() const;

	/**
	 * \brief Get the (axis-aligned) simulation domain.
	 */
	inline CPU void save(std::string file_name);

private:
	CPU void set_AABB(vec3 min, vec3 max);

	// These vectors are implicit 3D vectors, that represent a voxel grid
	std::vector<int> _mat_grid; // a voxel grid with nebula material codes. 
	std::vector<int> _tag_grid; // .. with PE tags
	std::vector<real> _e_grid; // .. with dissociation energies
	//std::vector<real> _dz_grid;  // .. with dz's
	std::vector<int> _species_grid;  // .. with electron species 
	
	real _voxel_size; // voxel size in nm
	int _size_x; // simulation domain size in x direction in voxels
	int _size_y; // .. y direction
	int _size_z; // .. z direction
	int _min_save_height = 0; // voxels with z-index < _min_save_height will not be saved
	int _max_save_height; //  voxels with z-index > _max_save_height will not be saved
	vec3 _AABB_min       = { 0, 0, 0 };
	vec3 _AABB_max       = { 0, 0, 0 };
	real _max_extent     = 0;

	friend struct detail::voxels_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "voxels.inl"

#endif // __GEOMETRY_VOXELS_H_
