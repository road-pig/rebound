import rebound
import unittest
import math
import rebound.data
import warnings
    
    
class TestIntegratorBS(unittest.TestCase):
    def test_bs_outersolarsystem(self):
        for eps in [1e-5, 1e-7, 1e-9, 1e-11]:
            sim = rebound.Simulation()
            rebound.data.add_outer_solar_system(sim)
            sim.integrator = "bs"
            sim.ri_bs.eps_rel = eps
            sim.ri_bs.eps_abs = eps
            e0 = sim.calculate_energy()
            sim.integrate(1000)
            e1 = sim.calculate_energy()
            self.assertLess(math.fabs((e0-e1)/e1),5*eps)

    def test_bs_high_e(self):
        for eps in [1e-7, 1e-9, 1e-11]:
            sim = rebound.Simulation()
            sim.add(m=1)
            sim.add(m=1e-3,a=1,e=0.9)
            sim.add(m=1e-3,a=6,e=0.9,f=0.5,omega=1.6)
            sim.integrator = "bs"
            sim.ri_bs.eps_rel = eps
            sim.ri_bs.eps_abs = eps
            e0 = sim.calculate_energy()
            sim.integrate(1000)
            e1 = sim.calculate_energy()
            self.assertLess(math.fabs((e0-e1)/e1),60*eps)
    
    def test_bs_inout(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        eps = 1e-6
        sim.integrator = "bs"
        sim.ri_bs.eps_rel = eps
        sim.ri_bs.eps_abs = eps
        sim.save("sim0.bin")
        sim1 = rebound.Simulation("sim0.bin")
        sim.integrate(100)
        sim1.integrate(100)
        sim1.save("sim1.bin")
        sim2 = rebound.Simulation("sim1.bin")
        sim.integrate(200)
        sim1.integrate(200)
        sim2.integrate(200)
        self.assertEqual(sim.particles[1].x, sim1.particles[1].x)
        self.assertEqual(sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(sim.particles[2].vx, sim1.particles[2].vx)
        self.assertEqual(sim.particles[2].vx, sim2.particles[2].vx)
        self.assertEqual(sim.t, sim1.t)
        self.assertEqual(sim.t, sim2.t)


    def test_bs_archive(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1)
        sim.add(m=1e-3,a=2,e=0.1)
        sim.automateSimulationArchive("test.sa",interval=10, deletefile=True)
        sim.integrate(100, exact_finish_time=0)
        sim1 = rebound.SimulationArchive("test.sa")[-1]
        sim.integrate(200, exact_finish_time=0)
        sim1.integrate(200, exact_finish_time=0)
        self.assertEqual(sim.particles[1].x, sim1.particles[1].x)
        self.assertEqual(sim.particles[2].vx, sim1.particles[2].vx)
        self.assertEqual(sim.t, sim1.t)
    
    def test_bs_change_particle_N(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1)
        sim.N_active = sim.N
        for i in range(100):
            sim.add(a=1.4+i*0.01)
        sim.integrate(10)
        for i in range(100):
            sim.add(a=2.4+i*0.01)
        sim.integrate(20)
        for i in range(50):
            sim.remove(i*2+2)
        sim.integrate(20)
        self.assertEqual(sim.N, 150+2)
    
    def test_bs_collide(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1,r=1,x=-3)
        sim.add(m=1,r=1,x=3)
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.integrate(20)
        self.assertEqual(sim.N, 1)
        self.assertEqual(sim.particles[0].m, 2.)
        self.assertEqual(sim.particles[0].x, 0.)
        self.assertEqual(sim.particles[0].vx, 0.)


if __name__ == "__main__":
    unittest.main()
