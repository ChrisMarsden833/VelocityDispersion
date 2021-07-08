#include "galaxy.h"

Galaxy::Galaxy(bool input_use_variable_IMF, float magnitude, float logsigma){
    use_variable_IMF = input_use_variable_IMF; 
    if(use_variable_IMF) this->get_IMF_params(magnitude, logsigma);
}

void Galaxy::get_IMF_params(float magnitude, float logsigma)
{
    if(-22.5 > magnitude) {
        if(logsigma > 2.4) R0 = 8., R8 = 4., R1 = 3.5;
        else R0 = 7., R8 = 3., R1  = 3.;       
    }
    else if( -21.5 > magnitude > -22.5){
        if(logsigma > 2.3) R0 = 5., R8 = 4., R1 = 3.;
        else R0 = 5., R8 = 3., R1 = 2.5;
    }
    else if(magnitude > -21.5)     {
        if(logsigma > 2.2) R0 = 5., R8 = 4., R1 = 3.;  
        else R0 = 3., R8 = 2.5, R1 = 2.;
    }
}

void Galaxy::ConstructBulge(float bulge_mag, float bulge_beta, float bulge_re, float bulge_n){
    //memberBulge = new Bulge(bulge_mag, bulge_re, bulge_n, bulge_beta, use_variable_IMF, R0, R1);
}

void Galaxy::ConstructHalo(float input_rho, float input_rs, const char * input_string){
    //memberHalo = new Halo(input_rho, input_rs, input_string);
}

void Galaxy::ConstructDisk(float disk_mass, float disk_h, float disk_i){
    //memberDisk = new Disk(disk_mass, disk_h, disk_i, memberBulge, memberHalo);
}



