import bfsimulator
import atomchip_rectwires
import cProfile

def run_sim():
    b_f_sim = bfsimulator.BFieldSimulator()
    b_f_sim.set_chip(atomchip_rectwires.atomchip_v4)
    b_f_sim.calc_trap_height()
    b_f_sim.plot_z()
    sim_results = b_f_sim.find_trap_freq()
    print 'x_trap : %2.0f um \ny_trap : %2.0f um \nz_trap : %2.0f um'%(b_f_sim.x_trap*1e6, b_f_sim.y_trap*1e6, b_f_sim.z_trap*1e6)
    print 'f_long : %2.0f Hz \nf_trans : %2.0f Hz \nf_z : %2.0f Hz'%(sim_results['f_long'], sim_results['f_trans'], sim_results['f_z'])
    # b_f_sim.plot_xy()
    
cProfile.run('run_sim()')