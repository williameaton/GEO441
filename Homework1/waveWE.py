from displacement import displacement_form
from vel_stress import vel_stress_form

def create_wave(type, model, BC_left, BC_right, solver=None, plot="v", label="Wave"):
    # Factory function that creates the type of wave wanted based on the different form of the wave equation being used
    # E.g. displacement or velocity-stress
        if type=="displacement":
            w = displacement_form(model, BC_left, BC_right, label)
            return w
        elif type=="vel-stress":
            w =  vel_stress_form(model, BC_left, BC_right, solver, label, plot)
            return w
        else:
            raise ValueError("Wave: type must be 'displacement' or 'vel-stress' ")