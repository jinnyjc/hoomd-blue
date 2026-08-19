// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "hoomd/mpcd/RejectionVirtualParticleFiller.h"
#ifdef ENABLE_HIP
#include "hoomd/mpcd/RejectionVirtualParticleFillerGPU.h"
#endif // ENABLE_HIP
#include "hoomd/mpcd/ParallelPlateGeometry.h"
#include "hoomd/mpcd/SphereGeometry.h"

#include "hoomd/SnapshotSystemData.h"
#include "hoomd/test/upp11_config.h"

HOOMD_UP_MAIN()

using namespace hoomd;

//! Check the virtual particles filled outside a sphere
/*!
 * Only the cells the surface can reach under grid shifting are filled, so the number of particles
 * has no closed form. Where they land is checked instead: a particle cannot sit farther from the
 * sphere than a marked cell can.
 */
template<class F>
void sphere_rejection_fill_basic_test(std::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    std::shared_ptr<SnapshotSystemData<Scalar>> snap(new SnapshotSystemData<Scalar>());
    snap->global_box = std::make_shared<BoxDim>(20.0);
    snap->particle_data.type_mapping.push_back("A");
    snap->mpcd_data.resize(1);
    snap->mpcd_data.type_mapping.push_back("A");
    snap->mpcd_data.type_mapping.push_back("B");
    snap->mpcd_data.position[0] = vec3<Scalar>(1, -2, 3);
    snap->mpcd_data.velocity[0] = vec3<Scalar>(123, 456, 789);
    std::shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));

    auto pdata = sysdef->getMPCDParticleData();
    auto cl = std::make_shared<mpcd::CellList>(sysdef, 1.0, true);
    // we should have no virtual particle in the system at this point.
    UP_ASSERT_EQUAL(pdata->getNVirtual(), 0);

    // create a spherical confinement of radius 5.0
    const Scalar r = 5.0;
    auto sphere = std::make_shared<const mpcd::SphereGeometry>(r, true);
    std::shared_ptr<Variant> kT = std::make_shared<VariantConstant>(1.5);
    std::shared_ptr<mpcd::RejectionVirtualParticleFiller<mpcd::SphereGeometry>> filler
        = std::make_shared<F>(sysdef, "B", 2.0, kT, sphere, 1000, 0);
    filler->setCellList(cl);

    /*
     * Test basic filling up for this cell list
     */
    unsigned int Nfill_0(0);
    filler->fill(0);
    // classification runs once inside the first fill, and the count is held for the second fill
    const unsigned int num_fill_cells_0 = filler->getNumFillCells();
        {
        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_vel(pdata->getVelocities(),
                                   access_location::host,
                                   access_mode::read);
        ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::read);

        // ensure first particle did not get overwritten
        UP_ASSERT_CLOSE(h_pos.data[0].x, Scalar(1), tol_small);
        UP_ASSERT_CLOSE(h_pos.data[0].y, Scalar(-2), tol_small);
        UP_ASSERT_CLOSE(h_pos.data[0].z, Scalar(3), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].x, Scalar(123), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].y, Scalar(456), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].z, Scalar(789), tol_small);
        UP_ASSERT_EQUAL(h_tag.data[0], 0);

        // check if the particles have been placed outside the confinement
        unsigned int N_out(0);
        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            // tag should equal index on one rank with one filler
            UP_ASSERT_EQUAL(h_tag.data[i], i);
            // type should be set
            UP_ASSERT_EQUAL(__scalar_as_int(h_pos.data[i].w), 1);

            Scalar3 pos = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
            const Scalar r2 = dot(pos, pos);
            if (r2 > r * r)
                ++N_out;
            }
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        Nfill_0 = N_out;
        }

    /*
     * Fill the volume again, which should approximately double the number of virtual particles
     */
    filler->fill(1);
    const unsigned int num_fill_cells_1 = filler->getNumFillCells();
        {
        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_vel(pdata->getVelocities(),
                                   access_location::host,
                                   access_mode::read);
        ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::read);

        // check if the particles have been placed outside the confinement
        unsigned int N_out(0);
        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            // tag should equal index on one rank with one filler
            UP_ASSERT_EQUAL(h_tag.data[i], i);
            // type should be set
            UP_ASSERT_EQUAL(__scalar_as_int(h_pos.data[i].w), 1);

            Scalar3 pos = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
            const Scalar r2 = dot(pos, pos);
            if (r2 > r * r)
                ++N_out;
            }
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        UP_ASSERT_GREATER(N_out, Nfill_0);

        // the second fill must reuse the same set of cells rather than classifying again
        UP_ASSERT_EQUAL(num_fill_cells_0, num_fill_cells_1);
        }

    /*
     * Test the average properties of the virtual particles.
     */
    // initialize variables for storing avg data
    Scalar N_avg(0);
    Scalar N_shell(0);
    Scalar3 vel_avg_net = make_scalar3(0, 0, 0);
    Scalar T_avg(0);

    // classification samples the cell grown by the maximum grid shift (0.5 in this case) on both
    // sides, so a marked cell can sit that much farther out. Therefore, a particle can be a body
    // diagonal away (sqrt(3)) for unit cells
    const Scalar r_out = r + std::sqrt(Scalar(3.0)) + 0.5;

    // repeat filling 10000 times
    unsigned int num_samples(10000);
    for (unsigned int t = 0; t < num_samples; ++t)
        {
        pdata->removeVirtualParticles();
        filler->fill(2 + t);

        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_vel(pdata->getVelocities(),
                                   access_location::host,
                                   access_mode::read);

        // local variables
        unsigned int N_out(0);
        Scalar temp(0);
        Scalar3 vel_avg = make_scalar3(0, 0, 0);

        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            const Scalar3 pos = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
            const Scalar3 vel = make_scalar3(h_vel.data[i].x, h_vel.data[i].y, h_vel.data[i].z);
            const Scalar r2 = dot(pos, pos);
            if (r2 > r * r)
                {
                ++N_out;
                if (r2 < r_out * r_out)
                    ++N_shell;
                }
            temp += dot(vel, vel);
            vel_avg += vel;
            }

        temp /= (3 * N_out);
        vel_avg_net += vel_avg / N_out;
        // Check whether all virtual particles are outside the sphere
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        N_avg += N_out;
        T_avg += temp;
        }
    N_avg /= num_samples;
    N_shell /= num_samples;
    T_avg /= num_samples;
    vel_avg_net /= num_samples;

    // every particle lies within the shell the marked cells occupy, so this should equal N_avg
    UP_ASSERT_EQUAL(N_shell, N_avg);

    UP_ASSERT_SMALL(vel_avg_net.x, tol_small);
    UP_ASSERT_SMALL(vel_avg_net.y, tol_small);
    UP_ASSERT_SMALL(vel_avg_net.z, tol_small);
    UP_ASSERT_CLOSE(T_avg, 1.5, tol_small);
    }

//! Check the virtual particles filled outside parallel plates
/*!
 * The walls sit at +/- separation/2 = +/- 2.5, so they cut the outermost cells exactly in half.
 * Classification must find those cells and no others, and they cover the whole solid region, so
 * the mean number of particles is known:
 *   N_exp = density * Lxz^2 * (Ly - separation)
 *
 * The walls are half a cell from a boundary, which is also the largest grid shift, so the same
 * cells are marked whether grid shifting is on or off and none of the expected values change.
 *
 * An xy tilt shears each y layer along x. It leaves the walls, the extent of a cell along y, and
 * the volume alone, so the expected values are the same as in the orthorhombic box.
 */
template<class F>
void plates_rejection_fill_test(std::shared_ptr<ExecutionConfiguration> exec_conf, Scalar xy)
    {
    const Scalar Lxz = 4.0;
    const Scalar Ly = 6.0;
    const Scalar a = 1.0;
    const Scalar separation = 5.0;
    const Scalar density = 10.0;
    const Scalar kT_val = 1.5;

    // 4 x 4 cells in each wall layer, two layers
    const unsigned int num_fill_cells_exp = 32;
    const Scalar N_exp = density * Lxz * Lxz * (Ly - separation);

    std::shared_ptr<SnapshotSystemData<Scalar>> snap(new SnapshotSystemData<Scalar>());
    snap->global_box = std::make_shared<BoxDim>(Lxz, Ly, Lxz);
    snap->global_box->setTiltFactors(xy, Scalar(0.0), Scalar(0.0));
    snap->particle_data.type_mapping.push_back("A");
    snap->mpcd_data.resize(1);
    snap->mpcd_data.type_mapping.push_back("A");
    snap->mpcd_data.type_mapping.push_back("B");
    snap->mpcd_data.position[0] = vec3<Scalar>(1, -2, 1);
    snap->mpcd_data.velocity[0] = vec3<Scalar>(123, 456, 789);
    std::shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));

    auto pdata = sysdef->getMPCDParticleData();
    auto cl = std::make_shared<mpcd::CellList>(sysdef, a, true);
    UP_ASSERT_EQUAL(pdata->getNVirtual(), 0);

    auto plates = std::make_shared<const mpcd::ParallelPlateGeometry>(separation, 0.0, true);
    std::shared_ptr<Variant> kT = std::make_shared<VariantConstant>(kT_val);
    std::shared_ptr<mpcd::RejectionVirtualParticleFiller<mpcd::ParallelPlateGeometry>> filler
        = std::make_shared<F>(sysdef, "B", density, kT, plates, 1000, 0);
    filler->setCellList(cl);

    /*
     * Test basic filling up for this cell list
     */
    unsigned int Nfill_0(0);
    filler->fill(0);
    const unsigned int num_fill_cells_0 = filler->getNumFillCells();
    UP_ASSERT_EQUAL(num_fill_cells_0, num_fill_cells_exp);
        {
        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_vel(pdata->getVelocities(),
                                   access_location::host,
                                   access_mode::read);
        ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::read);

        // ensure first particle did not get overwritten
        UP_ASSERT_CLOSE(h_pos.data[0].x, Scalar(1), tol_small);
        UP_ASSERT_CLOSE(h_pos.data[0].y, Scalar(-2), tol_small);
        UP_ASSERT_CLOSE(h_pos.data[0].z, Scalar(1), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].x, Scalar(123), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].y, Scalar(456), tol_small);
        UP_ASSERT_CLOSE(h_vel.data[0].z, Scalar(789), tol_small);
        UP_ASSERT_EQUAL(h_tag.data[0], 0);

        unsigned int N_out(0);
        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            // tag should equal index on one rank with one filler
            UP_ASSERT_EQUAL(h_tag.data[i], i);
            // type should be set
            UP_ASSERT_EQUAL(__scalar_as_int(h_pos.data[i].w), 1);

            const Scalar y = h_pos.data[i].y;
            if (std::abs(y) >= Scalar(0.5) * separation && std::abs(y) < Scalar(0.5) * Ly)
                ++N_out;
            }
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        Nfill_0 = N_out;
        }

    /*
     * Fill the volume again, which should approximately double the number of virtual particles
     */
    filler->fill(1);
    const unsigned int num_fill_cells_1 = filler->getNumFillCells();
        {
        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::read);

        unsigned int N_out(0);
        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            UP_ASSERT_EQUAL(h_tag.data[i], i);
            UP_ASSERT_EQUAL(__scalar_as_int(h_pos.data[i].w), 1);

            const Scalar y = h_pos.data[i].y;
            if (std::abs(y) >= Scalar(0.5) * separation && std::abs(y) < Scalar(0.5) * Ly)
                ++N_out;
            }
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        UP_ASSERT_GREATER(N_out, Nfill_0);
        // the second fill must reuse the same set of cells rather than classifying again
        UP_ASSERT_EQUAL(num_fill_cells_0, num_fill_cells_1);
        }

    /*
     * Test the average properties of the virtual particles
     */
    Scalar N_avg(0);
    Scalar3 vel_avg_net = make_scalar3(0, 0, 0);
    Scalar T_avg(0);
    unsigned int num_samples(20000);
    for (unsigned int t = 0; t < num_samples; ++t)
        {
        pdata->removeVirtualParticles();
        filler->fill(2 + t);

        ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> h_vel(pdata->getVelocities(),
                                   access_location::host,
                                   access_mode::read);

        unsigned int N_out(0);
        Scalar temp(0);
        Scalar3 vel_avg = make_scalar3(0, 0, 0);
        for (unsigned int i = pdata->getN(); i < pdata->getN() + pdata->getNVirtual(); ++i)
            {
            const Scalar y = h_pos.data[i].y;
            const Scalar3 vel = make_scalar3(h_vel.data[i].x, h_vel.data[i].y, h_vel.data[i].z);

            if (std::abs(y) >= Scalar(0.5) * separation && std::abs(y) < Scalar(0.5) * Ly)
                ++N_out;

            temp += dot(vel, vel);
            vel_avg += vel;
            }

        temp /= (3 * N_out);
        vel_avg_net += vel_avg / N_out;
        UP_ASSERT_EQUAL(N_out, pdata->getNVirtual());
        N_avg += N_out;
        T_avg += temp;
        }
    N_avg /= num_samples;
    T_avg /= num_samples;
    vel_avg_net /= num_samples;

    // the count per fill is Poisson, so the standard error of the mean is sqrt(N_exp/num_samples);
    // quoted at 4 sigma, and UP_ASSERT_CLOSE takes a percentage
    const Scalar tol_N = Scalar(100.0) * 4 * std::sqrt(N_exp / num_samples) / N_exp;
    UP_ASSERT_CLOSE(N_avg, N_exp, tol_N);

    UP_ASSERT_SMALL(vel_avg_net.x, tol_small);
    UP_ASSERT_SMALL(vel_avg_net.y, tol_small);
    UP_ASSERT_SMALL(vel_avg_net.z, tol_small);
    UP_ASSERT_CLOSE(T_avg, kT_val, tol_small);
    }

UP_TEST(sphere_rejection_fill_basic)
    {
    sphere_rejection_fill_basic_test<mpcd::RejectionVirtualParticleFiller<mpcd::SphereGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU));
    }

UP_TEST(plates_rejection_fill)
    {
    plates_rejection_fill_test<mpcd::RejectionVirtualParticleFiller<mpcd::ParallelPlateGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU),
        Scalar(0.0));
    }

UP_TEST(plates_rejection_fill_tilted)
    {
    plates_rejection_fill_test<mpcd::RejectionVirtualParticleFiller<mpcd::ParallelPlateGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU),
        Scalar(0.5));
    }
#ifdef ENABLE_HIP
UP_TEST(sphere_rejection_fill_basic_gpu)
    {
    sphere_rejection_fill_basic_test<mpcd::RejectionVirtualParticleFillerGPU<mpcd::SphereGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::GPU));
    }
UP_TEST(plates_rejection_fill_gpu)
    {
    plates_rejection_fill_test<
        mpcd::RejectionVirtualParticleFillerGPU<mpcd::ParallelPlateGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::GPU),
        Scalar(0.0));
    }

UP_TEST(plates_rejection_fill_tilted_gpu)
    {
    plates_rejection_fill_test<
        mpcd::RejectionVirtualParticleFillerGPU<mpcd::ParallelPlateGeometry>>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::GPU),
        Scalar(0.5));
    }
#endif // ENABLE_HIP
