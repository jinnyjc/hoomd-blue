// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "hoomd/mpcd/TriangulatedGeometry.h"
#include "hoomd/mpcd/TriangulatedGeometryFiller.h"

#include "hoomd/SnapshotSystemData.h"
#include "hoomd/test/upp11_config.h"

HOOMD_UP_MAIN()

using namespace hoomd;

//! Make parallel plates out of triangles
/*!
 * The plates are at y = +/- separation/2 with the fluid between them. Each plate is the cross
 * section of the box at its height, so in a box with an xy tilt it is displaced along x together
 * with the box. Every triangle is wound so that its normal points out of the fluid: +y for the
 * upper plate and -y for the lower plate.
 */
std::shared_ptr<const mpcd::TriangulatedGeometry>
makePlates(std::shared_ptr<SystemDefinition> sysdef, Scalar separation)
    {
    const BoxDim box = sysdef->getParticleData()->getGlobalBox();
    const Scalar Ly = box.getL().y;
    const Scalar f_hi = (Scalar(0.5) * separation + Scalar(0.5) * Ly) / Ly;
    const Scalar f_lo = (-Scalar(0.5) * separation + Scalar(0.5) * Ly) / Ly;

    std::vector<Scalar3> vertices;
    for (const Scalar f_y : {f_hi, f_lo})
        {
        vertices.push_back(box.makeCoordinates(make_scalar3(0, f_y, 0)));
        vertices.push_back(box.makeCoordinates(make_scalar3(0, f_y, 1)));
        vertices.push_back(box.makeCoordinates(make_scalar3(1, f_y, 0)));
        vertices.push_back(box.makeCoordinates(make_scalar3(1, f_y, 1)));
        }

    const std::vector<uint3> triangles = {// upper plate, normal +y
                                          make_uint3(0, 1, 3),
                                          make_uint3(0, 3, 2),
                                          // lower plate, normal -y
                                          make_uint3(4, 7, 5),
                                          make_uint3(4, 6, 7)};

    return std::make_shared<const mpcd::TriangulatedGeometry>(
        sysdef,
        static_cast<unsigned int>(vertices.size()),
        vertices.data(),
        static_cast<unsigned int>(triangles.size()),
        triangles.data(),
        Scalar(0.0),
        true);
    }

//! Check that two points are the same
void checkPoint(const Scalar3& p, const Scalar3& ref)
    {
    UP_ASSERT_SMALL(p.x - ref.x, tol_small);
    UP_ASSERT_SMALL(p.y - ref.y, tol_small);
    UP_ASSERT_SMALL(p.z - ref.z, tol_small);
    }

//! Closest point on a triangle
UP_TEST(triangulated_closest_point)
    {
    const Scalar3 a = make_scalar3(0, 0, 0);
    const Scalar3 b = make_scalar3(2, 0, 0);
    const Scalar3 c = make_scalar3(0, 2, 0);

    // face: the point projects onto the triangle, from either side
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(0.5, 0.5, 1.0), a, b, c),
               make_scalar3(0.5, 0.5, 0));
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(0.5, 0.5, -3.0), a, b, c),
               make_scalar3(0.5, 0.5, 0));

    // vertices
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(-1, -1, 0.5), a, b, c), a);
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(3, -0.5, 0), a, b, c), b);
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(-0.5, 3, 0), a, b, c), c);

    // edges ab, ac, and bc
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(1, -1, 0.3), a, b, c),
               make_scalar3(1, 0, 0));
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(-1, 1, 0), a, b, c),
               make_scalar3(0, 1, 0));
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(2, 2, 0.7), a, b, c),
               make_scalar3(1, 1, 0));

    // a point on the triangle is its own closest point
    checkPoint(mpcd::closestPointOnTriangle(make_scalar3(0.25, 0.5, 0), a, b, c),
               make_scalar3(0.25, 0.5, 0));
    }

//! Check if a point is inside or outside
UP_TEST(triangulated_is_outside)
    {
    auto exec_conf = std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU);
    std::shared_ptr<SnapshotSystemData<Scalar>> snap(new SnapshotSystemData<Scalar>());
    snap->global_box = std::make_shared<BoxDim>(10.0, 20.0, 10.0);
    snap->particle_data.type_mapping.push_back("A");
    snap->mpcd_data.resize(1);
    snap->mpcd_data.type_mapping.push_back("A");
    std::shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));

    auto plates = makePlates(sysdef, 15.0);
    UP_ASSERT_EQUAL(plates->getNumTotalTriangles(), 4);

    ArrayHandle<ShortReal3> h_verts(plates->getVertices(),
                                    access_location::host,
                                    access_mode::read);
    ArrayHandle<uint3> h_tris(plates->getTriangles(), access_location::host, access_mode::read);
    const std::vector<unsigned int> tris = {0, 1, 2, 3};

    auto isOutside = [&](const Scalar3& p)
    {
        return mpcd::isOutside(p,
                               tris.data(),
                               static_cast<unsigned int>(tris.size()),
                               h_verts.data,
                               h_tris.data);
    };

    // between the plates is fluid
    UP_ASSERT(!isOutside(make_scalar3(0, 0, 0)));
    UP_ASSERT(!isOutside(make_scalar3(1, 7, -2)));
    UP_ASSERT(!isOutside(make_scalar3(-3, -7, 4)));

    // beyond either plate is solid
    UP_ASSERT(isOutside(make_scalar3(1, 9, -2)));
    UP_ASSERT(isOutside(make_scalar3(-3, -8, 4)));

    // a point exactly on the surface is treated as fluid
    UP_ASSERT(!isOutside(make_scalar3(0, 7.5, 0)));

    // the diagonal edge shared by the two triangles of a plate gives the same answer
    UP_ASSERT(isOutside(make_scalar3(1, 8, 1)));
    UP_ASSERT(!isOutside(make_scalar3(1, 7, 1)));
    }

//! Check the virtual particles filled outside parallel plates made of triangles
/*!
 * This is the parallel plate test of the rejection filler with the plates replaced by a mesh,
 * so the expected values are the same. The walls sit at +/- separation/2 = +/- 7.5, so they cut
 * the cells containing them exactly in half. Classification must find those cells and no others,
 * so the mean number of particles is known:
 *   N_exp = density * (1/2) * num_fill_cells * a^3, where a is the cell size.
 *
 * The walls sit half a cell from a cell boundary, which is also the largest grid shift, so the
 * same cells are marked whether grid shifting is on or off and none of the expected values change.
 *
 * An xy tilt shears each y layer along x. It leaves the walls, the extent of a cell along y, and
 * the volume alone, so the expected values are the same as in the orthorhombic box.
 */
template<class F>
void triangulated_plates_fill_test(std::shared_ptr<ExecutionConfiguration> exec_conf, Scalar xy)
    {
    const Scalar Lxz = 10.0;
    const Scalar Ly = 20.0;
    const Scalar a = 1.0;
    const Scalar separation = 15.0;
    const Scalar density = 10.0;
    const Scalar kT_val = 1.5;

    // 10 x 10 cells in each wall layer, two layers
    const unsigned int num_fill_cells_exp = 200;
    const Scalar N_exp = density * Scalar(0.5) * num_fill_cells_exp;

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

    auto plates = makePlates(sysdef, separation);
    std::shared_ptr<Variant> kT = std::make_shared<VariantConstant>(kT_val);
    std::shared_ptr<mpcd::TriangulatedGeometryFiller> filler
        = std::make_shared<F>(sysdef, "B", density, kT, plates, 0);
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
        // the number of filled cells should not change between fills
        UP_ASSERT_EQUAL(num_fill_cells_0, num_fill_cells_1);
        }

    /*
     * Test the average properties of the virtual particles
     */
    Scalar N_avg(0);
    Scalar3 vel_avg_net = make_scalar3(0, 0, 0);
    Scalar T_avg(0);
    unsigned int num_samples(10000);
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
    // quoted at 4 sigma
    const Scalar tol_N = 4 * std::sqrt(N_exp / num_samples) / N_exp;
    UP_ASSERT_CLOSE(N_avg, N_exp, tol_N);

    // the mean velocity of one fill fluctuates with standard deviation sqrt(kT/N), so the average
    // over the samples has standard error sqrt(kT/(N*num_samples)); quoted at 4 sigma
    const Scalar tol_v = 4 * std::sqrt(kT_val / (N_exp * num_samples));
    UP_ASSERT_SMALL(vel_avg_net.x, tol_v);
    UP_ASSERT_SMALL(vel_avg_net.y, tol_v);
    UP_ASSERT_SMALL(vel_avg_net.z, tol_v);
    UP_ASSERT_CLOSE(T_avg, kT_val, tol_small);
    }

UP_TEST(triangulated_plates_fill)
    {
    triangulated_plates_fill_test<mpcd::TriangulatedGeometryFiller>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU),
        Scalar(0.0));
    }

UP_TEST(triangulated_plates_fill_tilted)
    {
    triangulated_plates_fill_test<mpcd::TriangulatedGeometryFiller>(
        std::make_shared<ExecutionConfiguration>(ExecutionConfiguration::CPU),
        Scalar(0.5));
    }
