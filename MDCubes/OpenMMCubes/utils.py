# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.


# State Report import

from __future__ import absolute_import
from __future__ import print_function

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

import simtk.openmm as mm
import simtk.unit as unit
import math
import time

# State Report End import

import parmed
from floe.api.orion import in_orion

from cuberecord import OELargeFile


class MDData(object):
    """
    This class is used to handle the MDData recovered
    from the Parmed structure.The class is designed to
    track changes in the pointed Parmed structure

    Notes
    -----
    Exposed variables:
        structure : Parmed structure
        positions : If present system atom positions otherwise None
        topology : Parmed topology
        box : If present box vectors otherwise None
        parameters : Parmed force field parameters
        velocities : If present system atom velocities otherwise None

    Examples
    --------
        mdData = MDData(parmed_structure)
        pos = mdData.positions
        vel = mdData.velocities
    """

    def __init__(self, parmed_structure):
        """
        Initialization function

        Parameters
        ----------
        parmed_structure : Parmed Structure object
            the parmed structure object
        """

        self.__parmed_structure__ = parmed_structure

        # Check atom positions
        if not self.__parmed_structure__.positions:
            raise RuntimeError('Atom positions are not defined')

    def __getattr__(self, attrname):
        if attrname == "structure":
            return self.__parmed_structure__
        elif attrname == "topology":
            return self.__parmed_structure__.topology
        elif attrname == "positions":
            # The returned object is an openmm Quantity with units
            return self.__parmed_structure__.positions
        elif attrname == "velocities":
            if self.__parmed_structure__.velocities is None:
                return None
            else:
                # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
                return self.__parmed_structure__.velocities * unit.angstrom/unit.picosecond
        elif attrname == "box":
            # The returned object is an openmm Quantity with units
            return self.__parmed_structure__.box_vectors
        elif attrname == "parameters":
            return parmed.ParameterSet.from_structure(self.__parmed_structure__)
        else:
            raise AttributeError('The required attribute is not defined: {}'.format(attrname))


def upload(filename):

    file_id = filename

    if in_orion():
        file_id = OELargeFile.create(filename)

    return file_id


def download(file_id):

    filename = file_id

    if in_orion():
        filename = file_id.retrieve()
        # if delete:
        #     file_id.delete()

    return filename


class StateDataReporterName(object):
    """
    This class has been adapted From OpenMM 7.1.1 to print the system name that is in process

    StateDataReporter outputs information about a simulation, such as energy and temperature, to a file.

    To use it, create a StateDataReporter, then add it to the Simulation's list of reporters.  The set of
    data to write is configurable using boolean flags passed to the constructor.  By default the data is
    written in comma-separated-value (CSV) format, but you can specify a different separator to use.
    """

    def __init__(self, file, reportInterval, system_name=None, step=False,
                 time=False, potentialEnergy=False, kineticEnergy=False,
                 totalEnergy=False, temperature=False, volume=False, density=False,
                 progress=False, remainingTime=False, speed=False, elapsedTime=False,
                 separator=',', systemMass=None, totalSteps=None):
        """Create a StateDataReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reportInterval : int
            The interval (in time steps) at which to write frames
        system_name : string=None
            The string that specify the system name
        step : bool=False
            Whether to write the current step index to the file
        time : bool=False
            Whether to write the current time to the file
        potentialEnergy : bool=False
            Whether to write the potential energy to the file
        kineticEnergy : bool=False
            Whether to write the kinetic energy to the file
        totalEnergy : bool=False
            Whether to write the total energy to the file
        temperature : bool=False
            Whether to write the instantaneous temperature to the file
        volume : bool=False
            Whether to write the periodic box volume to the file
        density : bool=False
            Whether to write the system density to the file
        progress : bool=False
            Whether to write current progress (percent completion) to the file.
            If this is True, you must also specify totalSteps.
        remainingTime : bool=False
            Whether to write an estimate of the remaining clock time until
            completion to the file.  If this is True, you must also specify
            totalSteps.
        speed : bool=False
            Whether to write an estimate of the simulation speed in ns/day to
            the file
        elapsedTime : bool=False
            Whether to write the elapsed time of the simulation in seconds to
            the file.
        separator : string=','
            The separator to use between columns in the file
        systemMass : mass=None
            The total mass to use for the system when reporting density.  If
            this is None (the default), the system mass is computed by summing
            the masses of all particles.  This parameter is useful when the
            particle masses do not reflect their actual physical mass, such as
            when some particles have had their masses set to 0 to immobilize
            them.
        totalSteps : int=None
            The total number of steps that will be included in the simulation.
            This is required if either progress or remainingTime is set to True,
            and defines how many steps will indicate 100% completion.
        """
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if (progress or remainingTime) and totalSteps is None:
            raise ValueError('Reporting progress or remaining time requires total steps to be specified')
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if file.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(file, 'wb', 0))
            elif file.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(file, 'w', 0)
            else:
                self._out = open(file, 'w')
        else:
            self._out = file
        self._system_name = system_name
        self._step = step
        self._time = time
        self._potentialEnergy = potentialEnergy
        self._kineticEnergy = kineticEnergy
        self._totalEnergy = totalEnergy
        self._temperature = temperature
        self._volume = volume
        self._density = density
        self._progress = progress
        self._remainingTime = remainingTime
        self._speed = speed
        self._elapsedTime = elapsedTime
        self._separator = separator
        self._totalMass = systemMass
        self._totalSteps = totalSteps
        self._hasInitialized = False
        self._needsPositions = False
        self._needsVelocities = False
        self._needsForces = False
        self._needEnergy = potentialEnergy or kineticEnergy or totalEnergy or temperature

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, self._needsPositions, self._needsVelocities, self._needsForces, self._needEnergy)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            print('#"%s"' % ('"'+self._separator+'"').join(headers), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialClockTime = time.time()
            self._initialSimulationTime = state.getTime()
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        # Check for errors.
        self._checkForErrors(simulation, state)

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values.
        print(self._separator.join(str(v) for v in values), file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, simulation, state):
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        box = state.getPeriodicBoxVectors()
        volume = box[0][0]*box[1][1]*box[2][2]
        clockTime = time.time()
        if self._system_name:
            values.append('%-10s' % self._system_name)
        if self._progress:
            values.append('%.1f%%' % (100.0*simulation.currentStep/self._totalSteps))
        if self._step:
            values.append(simulation.currentStep)
        if self._time:
            values.append(state.getTime().value_in_unit(unit.picosecond))
        if self._potentialEnergy:
            values.append(state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._kineticEnergy:
            values.append(state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole))
        if self._totalEnergy:
            values.append((state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole))
        if self._temperature:
            values.append((2*state.getKineticEnergy()/(self._dof*unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin))
        if self._volume:
            values.append(volume.value_in_unit(unit.nanometer**3))
        if self._density:
            values.append((self._totalMass/volume).value_in_unit(unit.gram/unit.item/unit.milliliter))
        if self._speed:
            elapsedDays = (clockTime-self._initialClockTime)/86400.0
            elapsedNs = (state.getTime()-self._initialSimulationTime).value_in_unit(unit.nanosecond)
            if elapsedDays > 0.0:
                values.append('%.3g' % (elapsedNs/elapsedDays))
            else:
                values.append('--')
        if self._elapsedTime:
            values.append(time.time() - self._initialClockTime)
        if self._remainingTime:
            elapsedSeconds = clockTime-self._initialClockTime
            elapsedSteps = simulation.currentStep-self._initialSteps
            if elapsedSteps == 0:
                value = '--'
            else:
                estimatedTotalSeconds = (self._totalSteps-self._initialSteps)*elapsedSeconds/elapsedSteps
                remainingSeconds = int(estimatedTotalSeconds-elapsedSeconds)
                remainingDays = remainingSeconds//86400
                remainingSeconds -= remainingDays*86400
                remainingHours = remainingSeconds//3600
                remainingSeconds -= remainingHours*3600
                remainingMinutes = remainingSeconds//60
                remainingSeconds -= remainingMinutes*60
                if remainingDays > 0:
                    value = "%d:%d:%02d:%02d" % (remainingDays, remainingHours, remainingMinutes, remainingSeconds)
                elif remainingHours > 0:
                    value = "%d:%02d:%02d" % (remainingHours, remainingMinutes, remainingSeconds)
                elif remainingMinutes > 0:
                    value = "%d:%02d" % (remainingMinutes, remainingSeconds)
                else:
                    value = "0:%02d" % remainingSeconds
            values.append(value)
        return values

    def _initializeConstants(self, simulation):
        """Initialize a set of constants required for the reports

        Parameters
        - simulation (Simulation) The simulation to generate a report for
        """
        system = simulation.system
        if self._temperature:
            # Compute the number of degrees of freedom.
            dof = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*unit.dalton:
                    dof += 3
            dof -= system.getNumConstraints()
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                dof -= 3
            self._dof = dof
        if self._density:
            if self._totalMass is None:
                # Compute the total system mass.
                self._totalMass = 0*unit.dalton
                for i in range(system.getNumParticles()):
                    self._totalMass += system.getParticleMass(i)
            elif not unit.is_quantity(self._totalMass):
                self._totalMass = self._totalMass*unit.dalton

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._system_name:
            headers.append('System')
        if self._progress:
            headers.append('Progress (%)')
        if self._step:
            headers.append('Step')
        if self._time:
            headers.append('Time (ps)')
        if self._potentialEnergy:
            headers.append('Potential Energy (kJ/mole)')
        if self._kineticEnergy:
            headers.append('Kinetic Energy (kJ/mole)')
        if self._totalEnergy:
            headers.append('Total Energy (kJ/mole)')
        if self._temperature:
            headers.append('Temperature (K)')
        if self._volume:
            headers.append('Box Volume (nm^3)')
        if self._density:
            headers.append('Density (g/mL)')
        if self._speed:
            headers.append('Speed (ns/day)')
        if self._elapsedTime:
            headers.append('Elapsed Time (s)')
        if self._remainingTime:
            headers.append('Time Remaining')
        return headers

    def _checkForErrors(self, simulation, state):
        """Check for errors in the current state of the simulation

        Parameters
         - simulation (Simulation) The Simulation to generate a report for
         - state (State) The current state of the simulation
        """
        if self._needEnergy:
            energy = (state.getKineticEnergy()+state.getPotentialEnergy()).value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                raise ValueError('Energy is NaN')
            if math.isinf(energy):
                raise ValueError('Energy is infinite')

    def __del__(self):
        if self._openedFile:
            self._out.close()
