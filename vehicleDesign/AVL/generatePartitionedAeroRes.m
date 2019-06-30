% close all;
clear;
clc;
format compact;
fclose all;

% create sample design
Partdsgn1 = avlDesignGeometryClass;

%% modify design
% wing parameters
% Input file names
Partdsgn1.input_file_name        = 'partDsgn1.avl'; % File name for .avl file
Partdsgn1.run_file_name          = 'testRunFile.run';
Partdsgn1.exe_file_name          = 'exeFile';
% Output file names
Partdsgn1.result_file_name       = 'partDsgn1_results';
Partdsgn1.lookup_table_file_name = 'partDsgn1_lookupTables';

% Name for design in the input file
Partdsgn1.design_name            = 'partDsgn1'; % String at top of input file defining the name

Partdsgn1.reference_point = [1.5;0;0];

Partdsgn1.wing_chord = 1.25;
Partdsgn1.wing_AR = 8;
Partdsgn1.wing_sweep = 8;
Partdsgn1.wing_dihedral = 2;
Partdsgn1.wing_TR = 0.8;
Partdsgn1.wing_incidence_angle = 0;
Partdsgn1.wing_naca_airfoil = '2412';
Partdsgn1.wing_airfoil_ClLimits = [-1.8 1.8];

Partdsgn1.wing_Nspanwise = 20;
Partdsgn1.wing_Nchordwise = 5;

Partdsgn1.h_stab_LE = 0.6*Partdsgn1.wing_chord*Partdsgn1.wing_AR;
Partdsgn1.h_stab_chord = 0.6*Partdsgn1.wing_chord;
Partdsgn1.h_stab_AR = 6;
Partdsgn1.h_stab_sweep = 10;
Partdsgn1.h_stab_dihedral = 0;
Partdsgn1.h_stab_TR = 0.9;
Partdsgn1.h_stab_naca_airfoil = '0012';
Partdsgn1.h_stab_airfoil_ClLimits = [-1.8 1.8];

Partdsgn1.v_stab_LE = Partdsgn1.h_stab_LE;
Partdsgn1.v_stab_chord = Partdsgn1.h_stab_chord;
Partdsgn1.v_stab_AR = 0.6*Partdsgn1.h_stab_AR;
Partdsgn1.v_stab_sweep = 15;
Partdsgn1.v_stab_TR = 0.8;
Partdsgn1.v_stab_naca_airfoil = '0012';
Partdsgn1.v_stab_airfoil_ClLimits = [-1.8 1.8];

%% run AVL
alpha_range = [-50 50];
n_steps = 51;

avlPartitioned(Partdsgn1,alpha_range,n_steps)

%% plot polars
plotPartitionedPolars(Partdsgn1.lookup_table_file_name)


