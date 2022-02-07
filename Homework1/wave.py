from displacement import displacement_form
from vel_stress import vel_stress_form

def create_wave(type, model, BC_left, BC_right, plot="v", label="Wave", IC=None):
        if type=="displacement":
            w = displacement_form(model, BC_left, BC_right, label)
            return w
        elif type=="vel-stress":
            w =  vel_stress_form(model, BC_left, BC_right, plot, label)
            return w
        else:
            raise ValueError("Wave: type must be 'displacement' or 'vel-stress' ")