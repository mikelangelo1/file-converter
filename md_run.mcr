# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Running an accurate molecular dynamics simulation in water with slow, normal or fast speed
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro sets up and runs a simulation. It can also continue a simulation that got interrupted.

# Parameter section - adjust as needed, but NOTE that some changes only take
# effect if you start an entirely new simulation, not if you continue an existing one. 
# ====================================================================================

# The structure to simulate must be present with a .pdb or .sce extension.
# If a .sce (=YASARA scene) file is present, the cell must have been added.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget 'c:\MyProject\1crn'

# pH at which the simulation should be run, by default physiological pH 7.4.
# To simulate in vacuo ('gas phase'), use ph='None' and pressurectrl='Off' further below,
# this will set functional groups to their neutral state found in vacuo.
# This setting affects the initially predicted protonation states only. If you want to compare
# simulations at different pH, make sure that you really have different protonation states.
ph=7.4

# The ion concentration as a mass fraction, here we use 0.9% NaCl (physiological solution)
ions='Na,Cl,0.9'

# Simulation temperature, which also serves as the random number seed (see Temp command).
# If you increase the temperature significantly by X%, you also need to reduce the timestep by X%
# by changing the 'tslist' that matches your speed below. If you run at a temperature that differs
# from 298K, you need to adapt the density below, look in the PressureCtrl documentation for a
# table of water densities as a function of temperature. Alternatively, change the pressurectrl
# below to activate Manometer1D mode.
temperature='298K'

# Water density in [g/ml], should match the temperature set above. If you do not know the proper
# density, make sure to enable 'Manometer1D' pressure control below.
density=0.997

# Pressure control mode
# Default: Rescale the cell such that residues named HOH reach the density specified above.
# This mode only makes sense if the solute is fully embedded in solvent, not for crystals or
# membranes. If your solvent is not water, create a single solvent molecule, set the property
# value of all atoms to the solvent density (Edit > Number > Property value), save it as
# YourStructure_solvent.yob, and enable the Manometer1D pressure control below. If your
# solvent is a mixture of several molecules please check the docs of the FillCellObj command.
pressurectrl='SolventProbe,Name=HOH,Density=(density)'

# Alternative: Uncomment below to calculate the pressure from the virial and
# uniformly rescale the cell to reach a pressure of 1 bar. Use this method if you
# do not know the correct density and your solute is still fully embedded in solvent.
#pressurectrl='Manometer1D,Pressure=1' 

# Alternative: Uncomment below to calculate the pressure from the virial and
# rescale the cell independently along each axis to reach a pressure of 1 bar.
# Use this method if the solute spans the entire cell (protein crystals...).
# See the PressureCtrl docs for other options, e.g. for membranes.
#pressurectrl='Manometer3D,Pressure=1' 

# Alternative: Do not control pressure, use NVT ensemble. Also use this for simulations in vacuo.
#pressurectrl='off'

# The format used to save the trajectories: YASARA 'sim', GROMACS 'xtc' or AMBER 'mdcrd'.
# If you don't pick 'sim', a single *.sim restart file will be saved too, since the other
# two formats don't contain velocities, only positions.
format='sim'

# Extension of the cell on each side around the solute in [A]
# '10' means that the cell will be 20 A larger than the protein.
# Cell settings only apply if you do not provide your own cell in a *.sce file.
extension=10

# Shape of the simulation cell: 'Cube', 'Cuboid' or 'Dodecahedron'.
# For long simulations that allow the solute to rotate freely, a dodecahedral cell
# is the fastest (watch the help movie 3.6 for details). For short simulations
# of elongated, non-spherical solutes, a rectangular 'Cuboid' box is likely faster,
# but it is then your responsibility to ensure that the simulation is really short enough
# and the solute does not rotate and interact with itself through periodic boundaries.
# Note that a dodecahedral cell needs much more memory than a cuboid cell, especially
# on GPUs. If you choose 'Dodecahedron' and plan to convert to XTC format later,
# click Options > Coordinate system > right-handed first. Note that Poisson-Boltzmann
# calculations (MM/PBSA) don't work in dodecahedral cells.
cellshape='Cube'

# The simulation speed, either 'slow' (2*1 fs timestep), 'normal' (2*1.25 fs timestep) or
# 'fast' (maximize performance with 2*2.5 fs timestep and constraints)
# Do not use 'fast' if you simulate incorrect molecules (that would not be stable in reality) 
# 'if !count speed' simply checks if variable 'speed' as been defined previously (e.g. by an including macro)
if !count speed
  speed='normal'

# Duration of the simulation, alternatively use e.g. duration=5000 to simulate for 5000 picoseconds
if !count duration
  duration='forever'

# The save interval for snapshots. Normally you don't need more than 500-1000 snapshots
# of your simulation, since that's the resolution limit of a typical figure in a journal.
if speed=='fast'
  # Fast speed, save simulation snapshots every 250000 fs, i.e. 250 ps.
  saveinterval=250000
else  
  # Slow or normal speed, save simulation snapshots every 100000 fs, i.e. 100 ps.
  saveinterval=100000

# Forcefield to use (these are now all YASARA commands, so no '=' used)
ForceField AMBER14

# Cutoff
Cutoff 8

# Cell boundary
Boundary periodic

# Use longrange coulomb forces (particle-mesh Ewald)
Longrange Coulomb

# Keep the solute from diffusing around and crossing periodic boundaries. Disable that for simulations of crystals.
CorrectDrift On

# Change the random seed to see how much your results change, or click Options > Random seed
#RandomSeed 1234567

# Normally no change required below this point
# ============================================

RequireVersion 15.1.1

# Treat all simulation warnings as errors that stop the macro
WarnIsError On

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

# When run as a macro in text mode, add configuration details to log file
if runWithMacro and ConsoleMode
  Processors

Clear
Console off
# Do we already have a scene with water or other solvent?
waterscene = FileSize (MacroTarget)_water.sce
solventscene = FileSize (MacroTarget)_solvent.sce

if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  # No scene with solvent present yet
  # Do we have a scene at all?
  scene = FileSize (MacroTarget).sce
  if scene
    LoadSce (MacroTarget)
    # Verify that the cell is present
    simcell = CountObj SimCell
    if !simcell
      RaiseError 'If you provide a scene, it must contain a simulation cell, but none was found in (MacroTarget).sce'
  else
    # No scene present, assume it's a PDB or YOB file
    for type in 'yob','pdb'
      size = FileSize (MacroTarget).(type)
      if size
        break
    if !size
      RaiseError 'Initial structure not found, expected (MacroTarget).pdb or .yob. Make sure to create a project directory and place the structure there'
    # Load structure
    Load(type) (MacroTarget)
    # In case user accidentally provided a YOb file with selected atoms
    Unselect
    # Align object with major axes to minimize cell size
    NiceOriAll
    # Delete long peptide bonds that bridge gaps in the structure, which tells CleanAll to add ACE/NME
    # capping groups (the structure of the missing residues could also be predicted, see LoadPDB docs).
    DelBond N,C,LenMin=5
    # Delete waters that are not involved in metal binding, to help the calculation of binding energies 
    DelRes Water with 0 arrows to all
    # Prepare the structure for simulation at the chosen pH
    CleanAll
    pH (ph)
    if Structure
      # Optimize the hydrogen-bonding network (more stable trajectories)
      OptHydAll
    # Create the simulation cell
    Cell Auto,Extension=(extension),Shape=(cellshape)
    SaveSce (MacroTarget)
  bnd = Boundary
  if bnd=='Wall'
    # The user supplied a cell with wall boundaries, we cannot use Experiment Neutralization
    ShowMessage "The simulation cell you created has wall boundaries, which reduces the simulation accuracy due to boundary effects..."
    Wait ContinueButton
    ShowMessage "You can click 'Simulation > Cell boundaries > Periodic' now to correct the problem, or 'Continue' immediately..."
    Wait ContinueButton
  if ph!='None'
    # Add water, user may have changed boundaries above
    bnd = Boundary
    if bnd=='Wall'
      # User really wants wall boundaries
      FillCellWater
    else
      Experiment Neutralization
        WaterDensity (density)
        pH (ph)
        Ions (ions)
        pKaFile (MacroTarget).pka
        Speed Fast
      Experiment On
      Wait ExpEnd
  # Do we have a solvent molecule with density stored as its property value?
  filename = '(MacroTarget)_solvent.yob'
  solventfound = FileSize (filename)
  if solventfound
    # Get rid of water molecules again, only need the counter ions to neutralize cell
    obj2 = ListObj Water
    DelRes Water
    # Load the solvent molecule and verify that the user set its density
    obj1 = LoadYOb (MacroTarget)_solvent
    dens = PropObj (obj1)
    if !dens
      RaiseError 'Please load the solvent molecule (filename), click Edit > Number > Property value > Obj X, choose the solvent density, then save the file again'
    CleanObj (obj1)
    # Fill the cell with solvent molecules
    FillCellObj (obj1),Density=(dens),BumpSum=4,RandomOri=100%
    # Join the solvent box to the counter ions and rename the object to 'Solvent'
    if obj1!=obj2
      JoinObj (obj)
    NameObj (obj2),Solvent
    SaveSce (MacroTarget)_solvent
  else
    # Save scene with water
    SaveSce (MacroTarget)_water
# Don't keep selected atoms, LoadXTC/LoadMDCRD would load only selected ones
Unselect

# Choose timestep and activate constraints
if speed=='fast'
  # Fast simulation speed
  # Constrain bonds to hydrogens
  FixBond all,Element H
  # Constrain certain bond angles involving hydrogens
  FixHydAngle all
  # Choose a multiple timestep of 2*2.5 = 5 fs
  # For structures with severe errors, 2*2 = 4 fs is safer (tslist=2,2)
  tslist=2,2.5
else
  # Slow or normal simulation speed
  # Remove any constraints
  FreeBond all,all
  FreeAngle all,all,all
  if speed=='slow'
    # Choose a multiple timestep of 2*1.00 = 2.0 fs
    tslist=2,1.0
  else
    # Choose a multiple timestep of 2*1.25 = 2.5 fs
    tslist=2,1.25
    # With this timestep, atoms may get too fast in very rare circumstances (only
    # in a specific protein, only once every few nanoseconds). The command below
    # slows down atoms moving faster than 13000 m/s. Such a 'random collision' every
    # few nanoseconds has no more impact than the random number seed. You can comment
    # it out for most proteins, or use the smaller timestep with speed 'slow' above:
    Brake 13000
# Update the pairlist every 10 (CPU) or 25 (GPU) steps
_,_,gpu = Processors
if gpu
  SimSteps Screen=25,Pairlist=25
else    
  SimSteps Screen=10,Pairlist=10
# Calculate total timestep, we want a float, so tslist2 is on the left side
ts=tslist2*tslist1
# Snapshots are saved every 'savesteps'
savesteps=saveinterval/ts
# Set final simulation parameters
TimeStep (tslist)
Temp (temperature)
# Check if user accidentally fixed some atoms
fixedlist() = ListAtom fixed
if count fixedlist
  if ConsoleMode
    FreeAll
  else
    MarkAtom (fixedlist1)
    ShowMessage '(count fixedlist) atoms are currently fixed. This will yield unrealistic trajectories, normally distances should be restrained instead, see user manual at Essentials > The 10 magic words > Bond. Click Simulation > Free > All if you agree...' 
    Wait ContinueButton

# Here you can make changes to the force field, just before the energy minimization starts

# Uncomment to fix certain atoms in space
#FixAtom Backbone Mol B

# Uncomment to add distance constraints
#AddSpring O Res Lys 80,H Res Glu 84,Len=1.9

# Uncomment to modify charges, e.g. let Trp 12 in Mol A lose an electron:
#ChargeRes Trp 12 Mol A,+1

# Alread a snapshot/trajectory present?
i=00000
if format=='sim'
  trajectfilename='(MacroTarget)(i).sim'
else  
  restartfilename='(MacroTarget).sim'
  trajectfilename='(MacroTarget).(format)'
  # Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
  old = FileSize (MacroTarget)(i).xtc
  if old
    RenameFile (MacroTarget)(i).xtc,(trajectfilename)
running = FileSize (trajectfilename)
if not running
  # Perform energy minimization
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  # And now start the real simulation
  Sim On
else
  # Simulation has been running before
  ShowMessage "Simulation has been running before, loading last snapshot..."
  if format=='sim'
    # Find and load the last SIM snapshot
    do
      i=i+1
      found = FileSize (MacroTarget)(i).sim
    while found
    i=i-1
    LoadSim (MacroTarget)(i)
    # Adjust savesteps to save snapshots in the same interval as previously
    if i>0
      t = Time
      savesteps=0+t/(ts*i)
  else
    # Do we have a restart file with atom velocities?
    found = FileSize (restartfilename)
    if found
      # Yes. First determine the savesteps if possible by loading the 2nd XTC/MDCrd snapshot
      last,t = Load(format) (trajectfilename),1
      if !last
        last,t = Load(format) (trajectfilename),2
        savesteps=0+t/ts
      # Then load the restart file
      LoadSim (restartfilename)
    else
      # No restart file found, load the last snapshot in the XTC/MDCrd trajectory
      do
        i=i+1
        last,t = Load(format) (trajectfilename),(i)
        ShowMessage 'Searching (format) trajectory for last snapshot, showing snapshot (i) at (0+t) fs'
        Sim Pause
        Wait 1
      while !last
      savesteps=0+t/(ts*(i-1))
      Sim Continue
HideMessage
  
# Set temperature and pressure control
TempCtrl Rescale
PressureCtrl (pressurectrl)

# And finally, make sure that future snapshots are saved
Save(format) (trajectfilename),(savesteps)
if format!='sim'
  # We additionally save a single SIM restart file with velocities
  SaveSim (restartfilename),(savesteps),Number=no

if duration=='forever'
  Console On
  if ConsoleMode
    # In the console, we need to wait forever to avoid a prompt for user input
    Wait forever
else
  measurements=0
  # Wait for given number of picoseconds, set 'duration' at the beginning.
  do
    # Tabulate properties you want to monitor during the simulation,
    # e.g. the speeds and velocity vectors of atoms 4, 5 and 7:
    # Tabulate SpeedAtom 4 5 7
    #
    # Or apply an acceleration (then you need to change Wait 10 to Wait 1 below)
    # AccelRes Protein,Y=1000
    #
    # Note that you can only read properties after each pairlist update.
    if ConsoleMode
      # In console model, display the ns/day
      Console On
    # For maximum frequency, use 'Wait 1' and - if needed - additionally reduce
    # the numbers at 'SimSteps' above.
    Wait 10
    Console Off
    measurements=measurements+1
    t = Time
    # Save PDB file of object 1 with picoseconds in the filename
    #SavePDB 1,(MacroTarget)_((0+t)//1000)
  while t<1000.*duration+1
  # Did we create a table with measurements?
  vallist() = Tab Default
  if count vallist
    # Yes, save the table
    SaveTab default,(MacroTarget)_duringsim,Format=Text,Columns=(count vallist/measurements),Header='Insert your own header here'
  Sim Off
# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit

  