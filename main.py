import numpy as np
import pylab
import copy
from qiskit import BasicAer
from qiskit.aqua import aqua_globals, QuantumInstance
from qiskit.aqua.algorithms import NumPyMinimunEigensolver, VQE
from qiskit.aqua.components.optimizers import SLSQP #ansatz
from qiskit.chemistry.components.initial_states import HartreeFock
from qisikt.chemistry.components.variational_forms import UCCSD
from qiskit.chemistry.drivers import PySCFDriver
from qiskit.chemistry.core import Hamiltonian, QubitMappingType 


Molecule = 'H .O .O -{O}; Li .O .O {O}'
Distances = np.arange(0.5,4.25,0.25)
vqe_energies = [] #ground satate energy calculated by vqe
hf_energies = []
exact_energies = [] #non optimized data

for i,d in enumerate(distances):
  print(step,i)

  #set up experiment
  driver = PySCFDriver(molecule.format(d/2), basis='sto3g')
  qmolecule = driver.run() #set up quantum molecule
  operator = Hamiltonian(qubit_mapping= QubitMappingType.PARITY, #hamiltonian represents the energy in quantum system
  two_qubit_reduction= True, freeze_core= True, orbital_reduction=[-3,-2])

  qubit_op, aux_ops = operator.run (qmolecule)

  #exact classical results
  exact_result = NumPyMinimumEigensolver (qubit_op,aux_operators = aux_ops).run()
  exact_result = operator.process_algorithm_result(exact_result)

  #VQE
  optimizier = SLSQP(maxiter=1000)
  initial_state = HartreeFock(operator.molecule_info['num_orbitals'], operator.molrcule_info['num_particals'], qubit_mapping=operator._qubit_mapping,two_qubit_reduction=operator._two_qubit_reduction)

  var_form = UCCSD(num_orbitals=operator.molecule_info['num_orbitals'], num_particles=operator.molecule_info['num_particles'], inisial_state=initial_state, qubit_mapping=operator._qubit_mapping, two_qubit_reduction=operator._two_qubit_reduction)

algo = VQE(qubit_op, var_form, optimizer, aux_operators=aux_ops)
vqe_result = algo.run(QuantumInstance(BasicAer.get_backend('statevector_simulator')))
vqe_result = operator.process_algorithm_result(vqe_result)

exact_energies.append(exact_result.energy)
vqe_energies.append(vqe_result.energy)
hf_energies.append(vqe_energies.hartree_fock_energy)

pylab.plot(distance, hf_energies, labal='Hartree_Fock')
pylab.plot(distances,vqe_energies, 'o', label='VQE')
pylab.plot(distances, exact_energies, 'x', label='Exact')

pylab.xlabel('Interatomic distance')
pylab.ylabel('Energies')
pylab.title('LiH Ground State Energy')
pylab1.legend(loc='upper right')






