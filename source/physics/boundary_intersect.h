#ifndef __BOUNDARY_INTERSECT_H_
#define __BOUNDARY_INTERSECT_H_

#include "../core/triangle.h"
#include "../geometry/voxels.h"

/**
 * \brief Material boundary intersection.
 *
 * See  (Kieft paper)
 *
 * \tparam quantum_transmission Consider probabilistic transmission
 * \tparam interface_refraction Consider refraction at interface
 * \tparam interface_absorption Empirical interface absorption (Kieft et al.
 *                                  doi:10.1088/0022-3727/41/21/215310)
 */
template<
	bool quantum_transmission = true,
	bool interface_refraction = true,
	bool interface_absorption = false,
	bool deposition = true>
struct boundary_intersect
{
	/**
	 * \brief Print diagnostic info
	 */
	nbl::geometry::voxels<false>* geometry;
	
	static void print_info(std::ostream& stream)
	{
		stream << std::boolalpha <<
			" * Material boundary crossing model\n"
			"   Options:\n"
			"     - Quantum mechanical transmission: " << quantum_transmission << "\n"
			"     - Interface refraction: " << interface_refraction << "\n"
			"     - Empirical interface absorption: " << interface_absorption << "\n";
	}

	/**
	 * \brief Perform intersection event.
	 */
	template<typename particle_manager, typename material_manager, bool gpu_flag>
	PHYSICS void execute(material_manager& material_mgr,
	                     particle_manager& particle_mgr, typename particle_manager::particle_index_t particle_idx,
	                     nbl::util::random_generator<gpu_flag>& rng) const
	{
		using material_index_t = typename material_manager::material_index_t;
		
		// Get particle data from particle manager
		auto this_particle = particle_mgr[particle_idx];
		//const triangle this_triangle = *(particle_mgr.get_last_triangle(particle_idx));
		material_index_t material_idx_in = particle_mgr.get_material_index(particle_idx);

		// Extract data from triangle pointer (code from luc)
		uint64_t isect_id = reinterpret_cast<uint64_t>(particle_mgr.get_last_triangle(particle_idx)); 
		material_index_t material_idx_out = reinterpret_cast<int32_t*>(&isect_id)[0];
		int voxel_side = reinterpret_cast<int32_t*>(&isect_id)[1];

		// Get particle direction, normal of the intersected triangle
		auto normalised_dir = normalised(this_particle.dir);
		vec3 last_triangle_normal;
		
			// determine the normal of the voxel side using its number in range (1..6)
		switch (voxel_side) 
		{
		case 1:
			last_triangle_normal = { 1, 0, 0 };
			break;

		case 2:
			last_triangle_normal = { -1, 0, 0 };
			break;

		case 3:
			last_triangle_normal = { 0, 1, 0 };
			break;

		case 4:
			last_triangle_normal = { 0, -1, 0 };
			break;

		case 5:
			last_triangle_normal = { 0, 0, 1 };
			break;

		case 6:
			last_triangle_normal = { 0, 0, -1 };
			break;
			
		default:	
			throw std::runtime_error("Invalid voxel side number!!!");
		}

		// Get angle between direction of motion and triangle
		const real cos_theta = dot_product(last_triangle_normal, normalised_dir);


		// determine the material index in the following way:
		//   material_idx_in represents the current material
		//   material_idx_out represents the material when passing through the interface
		/*material_index_t material_idx_in, material_idx_out;
		if (cos_theta > 0)
		{
			material_idx_in = this_triangle.material_in;
			material_idx_out = this_triangle.material_out;
		}
		else
		{
			material_idx_in = this_triangle.material_out;
			material_idx_out = this_triangle.material_in;
		}*/


		// manage special cases for electron detection, electron mirrors and terminators.
		//   DETECTOR always detects.
		//   DETECTOR_LT/GE50 detect under certain circumstances.
		//     If not detected, they pass through as if no intersection event has taken place.
		//
		
		
		switch (material_idx_out) {
		case material_manager::DETECTOR:			
			particle_mgr.detect(particle_idx);
			return;
		case material_manager::DETECTOR_LT50:
			if (this_particle.kin_energy < 50)
				particle_mgr.detect(particle_idx);
			return;
		case material_manager::DETECTOR_GE50:
			if (this_particle.kin_energy >= 50)
				particle_mgr.detect(particle_idx);
			return;
		case material_manager::TERMINATOR:
			particle_mgr.terminate(particle_idx);
			return;
		case material_manager::MIRROR:
			this_particle.dir = normalised_dir - 2*last_triangle_normal*cos_theta;
			particle_mgr[particle_idx] = this_particle;
			return;
		case material_manager::NOP:
			return;
		default:
			break;
		}
		

		// determine the change in energy `dU` (in eV) when passing through the interface
		// see thesis T.V. Eq. 3.136
		real dU = 0;
		if (material_mgr.is_physical(material_idx_in)) {
			dU -= material_mgr[material_idx_in].barrier;
		}
		if (material_mgr.is_physical(material_idx_out)) {
			dU += material_mgr[material_idx_out].barrier;
		}
		//std::clog << "\nParticle: " << particle_mgr.get_primary_tag(particle_idx) << "  " << voxel_side;
		// determine transmission probability (only if energy suffices)
		// see thesis T.V. Eq. 3.145
		if (this_particle.kin_energy*cos_theta*cos_theta + dU > 0)
		{
			const real s = sqrtr(1 + dU / (this_particle.kin_energy*cos_theta*cos_theta));
			const real T = (quantum_transmission
				? 4 * s / ((1 + s)*(1 + s))
				: 1);
			if (rng.unit() < T)
			{
				// if there is transmission, then adjust the kinetic energy,
				
				
				
				//this_particle.kin_energy += dU / 2;
				

				if (interface_refraction)
				{
					// determine the angle of refraction
					// see thesis T.V. Eq. 3.139
					this_particle.dir = (normalised_dir - last_triangle_normal * cos_theta)
						+ s * last_triangle_normal * cos_theta;
				}

				
				
				// update the current material index and EXIT.
								
				particle_mgr.set_material_index(particle_idx, material_idx_out);

				/*if (material_idx_out != material_manager::VACUUM)
				{
					this_particle.kin_energy += dU;
				}*/

				//this_particle.kin_energy += dU / 2;

				if (deposition)
				{
					// deposit a voxel of material 0
					
					
					if (material_idx_in == material_manager::VACUUM || material_idx_out == material_manager::VACUUM)
					{
						const real E = this_particle.kin_energy;
						const real dissociation_energy = 3.5; // The energy that an electron looses when dissociating an precursor molecule.
						real deposition_prob; // the probability of deposition according to a cross-section

						// Alman cross section 
						const real E_TH = 3.5; // dissosiation threshold energy in eV(waarden uit(2008))
						const real E_MAX = 18; // maximum dissosiation cross - section energy in eV
						const real LAMBDA_0 = 77; // lambda_0 in eV
						const real SIGMA_MAX = 1; // sigma_max

						// below are 4 dissosiation cross sections that can be switched on and off using multi line commenting

						/*if (this_particle.kin_energy <= E_TH) 
						{
							deposition_prob = 0;
						}
						else if (this_particle.kin_energy < E_MAX)
						{
							deposition_prob = SIGMA_MAX * (1 - ((E_MAX - this_particle.kin_energy) / std::pow(E_MAX - E_TH, 2)) );
						}
						else
						{
							deposition_prob = SIGMA_MAX * std::exp(-(this_particle.kin_energy - E_MAX) / LAMBDA_0);
						}*/

						// Alman-Winters cs
						/*if (this_particle.kin_energy <= E_TH)
						{
							deposition_prob = 0;
						}
						else if (this_particle.kin_energy < E_MAX)
						{
							deposition_prob = SIGMA_MAX * (1 - ((E_MAX - this_particle.kin_energy) / std::pow(E_MAX - E_TH, 2)) );
						}
						else
						{
							deposition_prob = E_MAX*std::log(this_particle.kin_energy / 6.621829941) / this_particle.kin_energy;
						}*/

						// Winters cross section
						/*const real E_MAX_w = 100;
						if(this_particle.kin_energy < 36.8)
						{
							deposition_prob = 0;
						}
						else
						{
							deposition_prob = 100*std::log(this_particle.kin_energy / 36.7879) / this_particle.kin_energy;
						}*/

						// Smith cross section
						/*if (E < 7)
						{
							deposition_prob = 0;
						}
						else if (E < 100)
						{
							deposition_prob = (1208 * (1 - 1 / E) + -1064 * std::pow((1 - 1 / E), 2) + -15.68 * std::log(E) + -797.4 * std::log(E) / E) / E;
						}
						else
						{
							deposition_prob = (-16540 * (1 - 1 / E) + 15970 * std::pow((1 - 1 / E), 2) + 108.8 * std::log(E) + 5885 * std::log(E) / E) / E;
						}*/

						deposition_prob = 0; // for sanity check only

						if (rng.unit() < deposition_prob)
						{
							/*vec3 dep_pos;
							if (material_idx_in == material_manager::VACUUM) // electron enters material from vacuum
							{
								dep_pos = -0.1 * last_triangle_normal + this_particle.pos; // deposition position
							}
							else if (material_idx_out == material_manager::VACUUM) // electron enters vacuum from material
							{
								dep_pos = 0.1 * last_triangle_normal + this_particle.pos; // deposition position
							}*/

							//std::clog <<  (this_particle.pos.x/0.3) << "   " << this_particle.pos.y/0.3 << "   " << this_particle.pos.z/0.3 << "   " << voxel_side;

							geometry->deposit(this_particle.pos, last_triangle_normal, 0, particle_mgr.get_primary_tag(particle_idx), this_particle.kin_energy, particle_mgr.get_species(particle_idx));

							this_particle.kin_energy -= dissociation_energy;
							//dU -= dissociation_energy;

							//particle_mgr.terminate(particle_idx); // After a deposition, the electron is terminated 

						}
					}
				}
				//this_particle.kin_energy += dU / 2;

				if (material_idx_out == material_manager::VACUUM)
				{
					particle_mgr.set_species(particle_idx, 3); // VE
					
				}
				this_particle.kin_energy += dU;
				
				
				particle_mgr[particle_idx] = this_particle;
				return;
			}
		}
		
		// surface absorption? (this is in accordance with Kieft & Bosch code)
		// note that the default behaviour has this feature disabled
		if (interface_absorption)
		{
			if (dU < 0 && rng.unit() < expr(1 + 0.5_r*this_particle.kin_energy / dU))
			{
				particle_mgr.terminate(particle_idx);
				return;
			}
		}

		// the only remaining case is total internal reflection
		this_particle.dir = normalised_dir - 2 * last_triangle_normal*cos_theta;
		particle_mgr[particle_idx] = this_particle;
	}
};

#endif // __BOUNDARY_INTERSECT_H_
