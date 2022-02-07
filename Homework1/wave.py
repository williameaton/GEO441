from displacement import displacement_form
from vel_stress import vel_stress_form
from initial_conditions import set_IC

def create_wave(type, model, BC_left, BC_right, plot="v", IC=set_IC):
        if type=="displacement":
            w = displacement_form(model, BC_left, BC_right, IC)
            return w
        elif type=="vel-stress":
            w =  vel_stress_form(model, BC_left, BC_right, IC, plot)
            return w
        else:
            raise ValueError("Wave: type must be 'displacement' or 'vel-stress' ")