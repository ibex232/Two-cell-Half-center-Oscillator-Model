from neuron import h
import matplotlib.pyplot as plt

class HCOCellTemplate(object):
  def __init__(self) -> None:
     self.create_soma()
     self.define_soma_dimensions()
     self.define_membrane_properties()
     self.setup_record()


### Create Cell Soma ###
  def create_soma(self):
    self.soma = h.Section(name= 'soma', cell=self)
    self.soma.nseg = 1

### Set Dimensions of Cell Soma ###
### Pi * 1000 * 9.99593 ~= 3.14e4 microns squared ###
  def define_soma_dimensions(self):
    self.soma.L = 1000 # 1mm in length
    self.soma.diam = 9.99593 # Same value as given in example file

  def define_membrane_properties(self):
    self.soma.cm = 1

    self.soma.insert('leak')
    self.soma.eleak = -60

    self.soma.insert('na')
    self.soma.ena = 50 # This is just taken from his rn,, check its right

    self.soma.insert('kdr')
    self.soma.ek = -80

    self.soma.insert('capool')
    self.soma.cao = 3
    self.soma.cai = 50e-6

    self.soma.insert('cas')
    self.soma.insert('ka')
    self.soma.insert('kca')
    self.soma.insert('cat')
    self.soma.insert('hyper')
    self.soma.eh = -20

    self.default_parameters = {
      'gbar_leak': 0.0004,
      'gbar_na': .12,    # (.1~.5)
      'gbar_kdr': .12,    # (.1~.5)
      'gbar_ka': 0.1,    # (.1~.5)
      'gbar_kca': 0.01,    # (.01~.05)
      'gbar_cas': 0.005,    # (.001~.01)
      'gbar_cat': 0.007,    # (.005~.01)
      'gbar_hyper': 0.0002,    # (.0001~.0003)
      'tauca_capool': 20,
      'fca_capool': 1.2
    }

    self.set_biophysics(**self.default_parameters)


#The following functions are taken and adapted from the provided HCO Cell template code
  def set_biophysics(self,**attributes):
    for param,value in attributes.items():
        if value is not None:
            setattr(self.soma,param,value)
        elif param in self.default_parameters.keys():
            setattr(self.soma,param,self.default_parameters[param])
    
  def get_biophysics(self,**attributes):
    for param in attributes.keys():
        attributes[param] = getattr(self.soma,param)
    return attributes
    
  def setup_record(self):
    self.t = h.Vector()
    self.t.record(h._ref_t)
    self.vars = ['ileak_leak','ina_na','ik_kdr',
                  'ica_cas','ica_cat','ik_ka','ik_kca',
                  'ih_hyper','v','cai']
    self.clrs = ['k','y','r','orange','brown','pink','g','c']
    self.record = {}
    for v in self.vars:
        vec = h.Vector()
        vec.record(getattr(self.soma(.5),'_ref_'+v))
        self.record[v] = vec
    
  def plot_vars(self,cellid=0,figsize=None):
    cellname = 'Cell B' if cellid>0 else 'Cell A'
    clr = 'r' if cellid>0 else 'b'
    t = self.t
    fig = plt.figure(figsize=figsize)
    axs = fig.subplots(3,1,sharex=True,gridspec_kw={'hspace':0.1})
    axs[0].set_title(cellname)
    axs[0].plot(t,self.record['v'],clr)
    axs[0].set_ylim(-90,60)
    axs[0].set_ylabel('Membrane Voltage (mV)')
    axs[2].plot(t,self.record['cai'],clr)
    axs[2].set_ylim(0,0.4)
    axs[2].set_ylabel('Calcium Pool (mM)')
    for i,v in enumerate(self.vars[:-2]):
        if getattr(self.soma,'gbar_'+v.split('_')[-1])>0:
            axs[1].plot(t,self.record[v],color=self.clrs[i],label=v)
    axs[1].legend(loc=1)
    axs[1].set_ylabel('Current (nA/cm$^2$)')
    axs[2].set_xlim(t[0],t[-1])
    axs[2].set_xlabel('Time (ms)')
    return fig, axs

