import os
from collections import defaultdict
import shutil

PERIODIC_TABLE = {
      1: "H",  2:"He",  3:"Li",  4:"Be",  5: "B",  6: "C",  7: "N",  8: "O",  9: "F", 10:"Ne",
     11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15: "P", 16: "S", 17:"Cl", 18:"Ar", 19: "K", 20:"Ca",
     21:"Sc", 22:"Ti", 23: "V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 30:"Zn",
     31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr", 37:"Rb", 38:"Sr", 39: "Y", 40:"Zr",
     41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn",
     51:"Sb", 52:"Te", 53: "I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La", 58:"Ce", 59:"Pr", 60:"Nd",
     61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy", 67:"Ho", 68:"Er", 69:"Tm", 70:"Yb",
     71:"Lu", 72:"Hf", 73:"Ta", 74: "W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg",
     81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At", 86:"Rn", 87:"Fr", 88:"Ra", 89:"Ac", 90:"Th",
     91:"Pa", 92: "U", 93:"Np", 94:"Pu", 95:"Am", 96:"Cm", 97:"Bk", 98:"Cf", 99:"Es",100:"Fm",
    101:"Md",102:"No",103:"Lr",104:"Rf",105:"Db",106:"Sg",107:"Bh",108:"Hs",109:"Mt",110:"Ds",
    111:"Rg",112:"Cn",113:"Nh",114:"Fl",115:"Mc",116:"Lv",117:"Ts",118:"Og"
}

START_TOKENS = {
    "xyz": ["Standard orientation:", "Input orientation:"],
    "freq": ["Harmonic frequencies (cm**-1)"],
    "energy": [" SCF Done:  "],
    "td": ["Excited states from <AA,BB:AA,BB> singles matrix:"],
    "pop": ["Population analysis using the SCF Density"],
    "spin": ["<Sx>="],
    "stable": ["Eigenvectors of the stability matrix:"],
}

END_TOKENS = {
    "xyz": ["Rotational constants (GHZ):", "Distance matrix (angstroms):"],
    "freq": ["Thermochemistry"],
    "energy": [" SCF Done: "],
    "td": ["SavETr:  write IOETrn="],
    "pop": ["Mulliken charges:"],
    "spin": ["<Sx>="],
    "stable": ["The wavefunction is"],
}

TARGET_TO_G16 = {
    "xyz": ["xyz"],
    "SCF_energy": ["energy"],
    "imaginary_frequencies": ["freq"],
    "<S**2>": ["spin"],
    "stable": ["stable"],
    "excitation_energy_s1": ["td"],
    "excitation_energy_t1": ["td"],
    "s1-t1_gap": ["td"],
}

def is_num(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True

def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    else:
        return True

def analyze_g16_output(id,base_dir,wfname,target_properties=["xyz"]):

    g16_properties = sorted(list(set([prop for target in target_properties for prop in TARGET_TO_G16[target]])))

    if not os.path.exists(base_dir / wfname / f"mol_{id}.log"):
        return False, dict()
    with open(base_dir / wfname / f"mol_{id}.log","r",encoding='utf-8') as f:
        lines = f.readlines()
    results = {k:"Fail" for k in g16_properties}
    reading_now = None
    reading_detail = None
    success = None
    reading = []
    reading_state = -1

    ln = 0 
    for line in lines:
        ln += 1
        # print(f"Line {ln}: {line.strip()}")
        # judge of success
        if "Normal termination of Gaussian 16" in line:
            success = True
        if "Error termination" in line:
            success = False

        # activate
        if reading_now is None:
            for prop in g16_properties:
                if any(token in line for token in START_TOKENS[prop]):
                    # print(f"Start reading {prop}. (Line {ln})")
                    reading_now = prop
                    break
        
        # capture
        if reading_now == "xyz":
            l = [w for w in line.replace("\n","").split(" ") if w != ""]
            if len(l) != 0 and is_num(l[0]):
                reading.append(tuple([PERIODIC_TABLE[int(l[1])]]+[f'{float(val): .6f}' for val in l][3:6]))

        elif reading_now == "energy":
            reading.append(float([w for w in line.split(" ") if w != ""][4]))

        elif reading_now == "freq":
            line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
            if len(line_contents) == 0:
                pass
            elif is_int(line_contents[0]) and len(line_contents) <= 3:
                reading_state = len(line_contents)
                for _ in range(reading_state):
                    reading.append(dict())
                    reading[-1]["displacement"] = []
            elif len(line_contents) == reading_state:
                for i in range(reading_state):
                    reading[-reading_state+i]["symmetry"] = line_contents[i]
            elif line_contents[0] == "Frequencies":
                for i in range(reading_state):
                    reading[-reading_state+i]["frequencies"] = float(line_contents[i+2])
            elif line_contents[0] == "Red.":
                for i in range(reading_state):
                    reading[-reading_state+i]["red_mass"] = float(line_contents[i+3])
            elif line_contents[0] == "Frc":
                for i in range(reading_state):
                    reading[-reading_state+i]["force_const"] = float(line_contents[i+3])
            elif line_contents[0] == "IR":
                for i in range(reading_state):
                    reading[-reading_state+i]["IR_int"] = float(line_contents[i+3])
            elif is_int(line_contents[0]):
                for i in range(reading_state):
                    reading[-reading_state+i]["displacement"].append(tuple([float(val) for val in line_contents[i*3+2:i*3+5]]))

        elif reading_now == "td":
            if reading_detail is None:
            #     reading_detail = "td_iter"

            # if reading_detail == "td_iter":
            #     if "converged" in line:
            #         pass
            #     elif "Root" in line:
            #         line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
            #         if len(reading) < int(line_contents[1]):
            #             for _ in range(int(line_contents[1]) - len(reading)):
            #                 reading.append(dict())
            #         reading[int(line_contents[1])-1]["energy"] = float(line_contents[3])
                if "Ground to excited state transition electric dipole moments (Au):" in line:
                    reading_detail = "transition_electric_dipole_moments"

            elif reading_detail == "transition_electric_dipole_moments":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    if len(reading) < int(line_contents[0]):
                        for _ in range(int(line_contents[0]) - len(reading)):
                            reading.append(dict())
                    reading[int(line_contents[0])-1]["transition_electric_dipole_moments"] = {
                        "x": float(line_contents[1]),
                        "y": float(line_contents[2]),
                        "z": float(line_contents[3]),
                        "Dipole_Strength": float(line_contents[4]),
                        "Oscillator_Strength": float(line_contents[5]),
                    }
                elif "Ground to excited state transition velocity dipole moments (Au):" in line:
                    reading_detail = "transition_velocity_dipole_moments"

            elif reading_detail == "transition_velocity_dipole_moments":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["transition_velocity_dipole_moments"] = {
                        "x": float(line_contents[1]),
                        "y": float(line_contents[2]),
                        "z": float(line_contents[3]),
                        "Dipole_Strength": float(line_contents[4]),
                        "Oscillator_Strength": float(line_contents[5]),
                    }
                elif "Ground to excited state transition magnetic dipole moments (Au):" in line:
                    reading_detail = "transition_magnetic_dipole_moments"

            elif reading_detail == "transition_magnetic_dipole_moments":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["transition_magnetic_dipole_moments"] = {
                        "x": float(line_contents[1]),
                        "y": float(line_contents[2]),
                        "z": float(line_contents[3]),
                    }
                elif "Ground to excited state transition velocity quadrupole moments (Au):" in line:
                    reading_detail = "transition_velocity_quadrupole_moments"

            elif reading_detail == "transition_velocity_quadrupole_moments":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["transition_velocity_quadrupole_moments"] = {
                        "xx": float(line_contents[1]),
                        "yy": float(line_contents[2]),
                        "zz": float(line_contents[3]),
                        "xy": float(line_contents[4]),
                        "xz": float(line_contents[5]),
                        "yz": float(line_contents[6]),
                    }
                elif "Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)" in line:
                    reading_detail = "rotatory_strengths"
            
            elif reading_detail == "rotatory_strengths":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["rotatory_strengths"] = {
                        "xx": float(line_contents[1]),
                        "yy": float(line_contents[2]),
                        "zz": float(line_contents[3]),
                        "R_velocity": float(line_contents[4]),
                        "E-M_angle": float(line_contents[5]),
                    }
                elif "Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)" in line:
                    reading_detail = "rotatory_strengths2"
            
            elif reading_detail == "rotatory_strengths2":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["rotatory_strengths2"] = {
                        "xx": float(line_contents[1]),
                        "yy": float(line_contents[2]),
                        "zz": float(line_contents[3]),
                        "R_length": float(line_contents[4]),
                    }
                elif "1/2[<0|del|b>*<b|r|0> + (<0|r|b>*<b|del|0>)*] (Au)" in line:
                    reading_detail = "rotatory_strengths3"
            
            elif reading_detail == "rotatory_strengths3":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if len(line_contents) == 0 or line_contents[0] == "state":
                    pass
                elif is_int(line_contents[0]):
                    reading[int(line_contents[0])-1]["rotatory_strengths3"] = {
                        "x": float(line_contents[1]),
                        "y": float(line_contents[2]),
                        "z": float(line_contents[3]),
                        "Dipole_Strength": float(line_contents[4]),
                        "Oscillator_Strength": float(line_contents[5]),
                    }
                elif "Excitation energies and oscillator strengths:" in line:
                    reading_detail = "td_final"
            
            elif reading_detail == "td_final":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if len(line_contents) == 0:
                    pass
                elif line_contents[0] == "Excited" and line_contents[1] == "State":
                    reading_state = int(line_contents[2][:-1])
                    reading[reading_state-1]["detail"] = {
                        "Symmetry": line_contents[3],
                        "Energy": float(line_contents[4]),
                        "Wavelength": float(line_contents[6]),
                        "Oscillator_Strength": float(line_contents[8][2:]),
                        "<S**2>": float(line_contents[9][7:]),
                    }
                    reading[reading_state-1]["configuration"] = []
                elif is_int(line_contents[0]):
                    reading[reading_state-1]["configuration"].append({
                        "transition": "".join(line_contents[:-1]),
                        "coefficient": float(line_contents[-1]),
                    })
                elif "Total Energy" in line:
                    if reading_state != -1:
                        reading[reading_state-1]["TD_total_energy"] = float(line_contents[4])
        
        elif reading_now == "pop":
            if reading_detail is None:
                if "Orbital symmetries:" in line:
                    reading_detail = "symmetries"
                elif "Alpha  occ. eigenvalues --" in line:
                    reading_detail = "eigenvalues"

            if reading_detail == "symmetries":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if len(line_contents) == 0:
                    pass
                elif line_contents[0] == "Occupied":
                    reading_state = "occupied"
                    for i in range(1, len(line_contents)):
                        reading.append(dict())
                        reading[-1]["symmetry"] = line_contents[i]
                        reading[-1]["occupied"] = True
                elif line_contents[0] == "Virtual":
                    reading_state = "virtual"
                    for i in range(1, len(line_contents)):
                        reading.append(dict())
                        reading[-1]["symmetry"] = line_contents[i]
                        reading[-1]["occupied"] = False
                elif line_contents[0][0] == "(":
                    for i in range(len(line_contents)):
                        reading.append(dict())
                        reading[-1]["symmetry"] = line_contents[i]
                        reading[-1]["occupied"] = (reading_state == "occupied")
                elif "The electronic state " in line:
                    reading_detail = "eigenvalues"
                    reading_state = 0
            
            elif reading_detail == "eigenvalues":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if line_contents[2] == "eigenvalues":
                    for i in range(4, len(line_contents)):
                        reading[reading_state]["eigenvalue"] = float(line_contents[i])
                        reading_state += 1
                elif "Molecular Orbital Coefficients:" in line:
                    reading_detail = "molecular_orbital_coefficients"
                elif "Condensed to atoms" in line:
                    reading_detail = "condensed_to_atoms"
            
            elif reading_detail == "molecular_orbital_coefficients":
                line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
                if all([is_int(w) for w in line_contents]):
                    reading_state = [int(w)-1 for w in line_contents]
                    for orb in reading_state:
                        reading[orb]["coefficients"] = []
                elif is_int(line_contents[0]):
                    for i, orb in enumerate(reading_state[::-1]):
                        reading[orb]["coefficients"].append(float(line_contents[-(i+1)]))
                elif "Density Matrix:" in line:
                    reading_detail = "density_matrix"
            
            elif reading_detail == "density_matrix":
                if "Full Mulliken population analysis:" in line:
                    reading_detail = "full_mulliken_population_analysis"

            elif reading_detail == "full_mulliken_population_analysis":
                if "Gross orbital populations:" in line:
                    reading_detail = "gross_orbital_populations"
            
            elif reading_detail == "gross_orbital_populations":
                if "Mulliken charges:" in line:
                    reading_detail = "mulliken_charges"

        elif reading_now == "spin":
            line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
            reading.append({
                "Sx": float(line_contents[1]),
                "Sy": float(line_contents[3]),
                "Sz": float(line_contents[5]),
                "S**2": float(line_contents[7]),
                "S": float(line_contents[9]),
            })
        
        elif reading_now == "stable":
            line_contents = [w for w in line.replace("\n","").split(" ") if w != ""]
            if len(line_contents) != 0 and line_contents[0] != "Eigenvectors":
                reading.append({"Eigenvalue": float(line_contents[4]), "<S**2>": float(line_contents[5][7:])})
        
        # deactivate
        if reading_now in END_TOKENS and any(token in line for token in END_TOKENS[reading_now]):
            # print(f"Finished reading {reading_now}. (Line {ln})")
            result = None
            if reading_now == "xyz":
                result = tuple(reading)
            elif reading_now == "freq":
                result = reading
            elif reading_now == "energy":
                result = reading[0]
            elif reading_now == "td":
                result = reading
            elif reading_now == "pop":
                result = reading
            # elif reading_now == "orbital":
            #     result = tuple(reading)
            # elif reading_now == "dipole":
            #     result = reading[0]
            elif reading_now == "spin":
                result = reading[0]
            elif reading_now == "stable":
                result = reading
            results[reading_now] = result
            reading_now = None
            reading = []
            reading_detail = None
            reading_state = -1

    return success,results

def output_to_property(results, target_properties):
    properties = dict()

    for target in target_properties:
        if target == "xyz":
            properties[target] = results["xyz"]

        elif target == "SCF_energy":
            properties[target] = results["energy"]

        elif target == "imaginary_frequencies":
            properties[target] = results["freq"][0]["frequencies"] > 0

        elif target == "<S**2>":
            properties[target] = results["spin"]["S**2"]

        elif target == "stable":
            properties[target] = results["stable"][0]["Eigenvalue"] > 0
        
        elif target == "excitation_energy_s1":
            for td_results in results["td"]:
                if "Singlet" in td_results["detail"]["Symmetry"]:
                    properties[target] = td_results["detail"]["Energy"]
                    break
        
        elif target == "excitation_energy_t1":
            for td_results in results["td"]:
                if "Triplet" in td_results["detail"]["Symmetry"]:
                    properties[target] = td_results["detail"]["Energy"]
                    break
        
        elif target == "s1-t1_gap":
            s1_energy = None
            t1_energy = None
            for td_results in results["td"]:
                if "Singlet" in td_results["detail"]["Symmetry"]:
                    if s1_energy is None:
                        s1_energy = td_results["detail"]["Energy"]
                elif "Triplet" in td_results["detail"]["Symmetry"]:
                    if t1_energy is None:
                        t1_energy = td_results["detail"]["Energy"]
            if s1_energy is not None and t1_energy is not None:
                properties[target] = s1_energy - t1_energy
            else:
                properties[target] = None

    return properties
