#This file will contain the functions to describe magnetospheric boundaries and functions to describe the emissivity. 

#You wouldn't run this file directly. 

import numpy as np 
import warnings
warnings.filterwarnings('ignore')

#SHUE BOUNDARY MODELS. 
######################

def shue_func(theta, phi, r0, ay, az):
        '''This is the 3D Shue model defined in Jorgensen et al. (2019)
        
        Parameters
        ----------
        theta (rad) and phi (rad)
        r0 - subsolar magnetopause distance
        ay - alpha y parameter
        az - alpha z parameter 

        Returns
        -------
        r - radial distance at the angles theta and phi 
        '''

        ry = r0*((2/(1+np.cos(theta)))**ay)
        rz = r0*((2/(1+np.cos(theta)))**az)

        r = (ry*rz)/(((rz*np.cos(phi))**2 + (ry*np.sin(phi))**2)**0.5)

        return r 

#def shue_adapted_func(theta, phi, rmp0, dr, ay, az):
#        '''This is an adapted 3D Shue model that doesn't use the bowshock subsolar distance but adds to the magnetopause distance instead. Used for CMEM2.
        
#        Parameters
#        ----------
#        theta (rad) and phi (rad) 
#        rmp0 - subsolar magnetopause distance (calculated with get_rmp0)
#        dr - delta r parameter. distance between magnetopause and bowshock along subsolar line. 
#        ay - alpha y parameter
#        az - alpha z parameter 
        
#        Returns
#        -------
#        r - radial distance at the angles theta and phi
#        '''
        
#        ry = (rmp0+dr)*((2/(1+np.cos(theta)))**ay)
#        rz = (rmp0+dr)*((2/(1+np.cos(theta)))**az)
        
#        r = (ry*rz)/(((rz*np.cos(phi))**2 + (ry*np.sin(phi))**2)**0.5)

#        return r 
  
def shue_func_simple(theta, rmp0, dr, p1, dp1):
    '''This is a simplified bowshock model like the original Shue model that is not axisymmetric.''' 
    
    return (rmp0+dr)*((2/(1+np.cos(theta)))**(p1+dp1)) 

def shue_func_simple_cmem2f(theta, rmp0, p1, dp1):
    '''This is a simplified bowshock model like the original Shue model that is not axisymmetric.''' 
    
    return (1.626*rmp0-3.166)*((2/(1+np.cos(theta)))**(p1+dp1)) 

#MAGNETOPAUSE BOUNDARY MODELS.
##############################
     
def lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0=1, p1=1, p2=1, p3=1):
        '''This function will work out r using the lin model. 
        
        Parameters
        ----------
        theta (rad) - Shue coords.
        phi (rad) - Shue coords. 
        a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
        #dipole - dipole tilt angle (rad)
        #pd - dynamic pressure in nPa
        #pm - magnetic pressure in nPa 
        #bz - IMF bz component in nT 
        p - parameter scaling factors. 
            p0 scales r0
            p1 scales flaring parameter beta 
            p2 scales indentation parameter Q (cusp depth) 
            p3 scales d in indentation shape (cusp shape/width)
            '''

        # Get coefficients if for some reason, they have not already been calculated. 
        #if r0_lin is None: 
        #    get_lin_coeffs(dipole, pd, pm, bz)
        
        # Get phi-n and phi-s.
        phi_n = np.arccos((np.cos(theta)*np.cos(theta_n)) + (np.sin(theta)*np.sin(theta_n)*np.cos(phi-(np.pi/2.))))
        phi_s = np.arccos((np.cos(theta)*np.cos(theta_s)) + (np.sin(theta)*np.sin(theta_s)*np.cos(phi-(3*np.pi/2.))))

        # Get f. 
        f = (np.cos(theta/2) + a[5]*np.sin(2*theta)*(1-np.exp(-theta)))**(p1*(beta_c[0] + beta_c[1]*np.cos(phi) + beta_c[2]*np.sin(phi) + beta_c[3]*(np.sin(phi)**2)))

        # Get Q. 
        Q = p2*c*np.exp(p3*dn*(phi_n**a[21])) + p2*c*np.exp(p3*ds*(phi_s**a[21]))

        # Get r. 
        r = p0*r0_lin*f + Q

        return r 

def lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, p0, p1, p2, p3):
        '''This will calculate the simplified version of the boundary model that does not contain any solar wind dependences.
        
        Parameters
        ----------
        theta (rad) - Shue coords.
        phi (rad) - Shue coords. 
        a, dn, ds, theta_n, theta_s, r0_lin - Simplified Lin coefficients in model. 
        p - parameter scaling factors. 
            p0 scales r0
            p1 scales flaring parameter beta 
            p2 scales indentation parameter Q (cusp depth) 
            p3 scales d in indentation shape (cusp shape/width)''' 
        
        # Get phi-n and phi-s.
        phi_n = np.arccos((np.cos(theta)*np.cos(theta_n)) + (np.sin(theta)*np.sin(theta_n)*np.cos(phi-(np.pi/2.))))
        phi_s = np.arccos((np.cos(theta)*np.cos(theta_s)) + (np.sin(theta)*np.sin(theta_s)*np.cos(phi-(3*np.pi/2.))))
        
        # Get f. 
        f = (np.cos(theta/2) + a[5]*np.sin(2*theta)*(1-np.exp(-theta)))**(-p1)

        # Get Q. 
        Q = p2*np.exp(p3*dn*(phi_n**a[21])) + p2*np.exp(p3*ds*(phi_s**a[21]))
        
        # Get r. 
        r = p0*f - Q
        
        return r 

        
#GET COEFFICIENTS FOR BOUNDARY MODELS. 
######################################


def get_rmp0(a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p2, p3):
        '''This will get the subsolar magnetopause distance for the current parameters. Original, hence not subtracting Q0.  
        
        Parameters
        ----------
        Lin parameters from get_lin_coeffs 
        p0
        p2
        p3
        
        '''
        
        Q0 = c*(np.exp(p3*dn*(theta_n**a[21])) + np.exp(p3*ds*(theta_s**a[21]))) 
        
        rmp0 = p0*r0_lin + p2*Q0 
        
        return rmp0 

def get_rmp0_simple(a, dn, ds, theta_n, theta_s, p0, p2, p3):
        '''This will get the subsolar magnetopause for the simplified CMEM model without solar wind dependence. 
        '''
        
        Q0 = p2*(np.exp(p3*dn*(theta_n)**a[21])) + p2*(np.exp(p3*dn*(theta_s)**a[21]))
        
        rmp0 = p0 - Q0
        
        return rmp0
        
        
#def get_k(maxIx, maxdIx, sigmoid='exp'):
#        '''This will create an estimate for the boundary parameter k based on values in the emissivity cube.
        
#        Parameters
#        ----------
#        maxIx - Radial position of maximum intensity along sun-earth line. 
#        maxdIx - Radial position of maximum gradient in intensity along sun-earth line. 
#        sigmoid - 'exp' or 'arctan'. Will return different values for each. 
        
#        ''' 

#        if sigmoid == 'exp':
#            k = 2*np.log(99)/abs(maxIx-maxdIx) 
#            return k     
#        elif sigmoid == 'arctan':
#            k = (2*np.tan(0.49*np.pi))/abs(maxIx-maxdIx) 
#            return k 
#        else:
#            raise ValueError("Invalid sigmoid function chosen: 'exp' or 'arctan'")
 
def get_lin_coeffs(dipole, pd, pm, bz):
        '''This gets the value of r0 in the Lin et al. (2010) model, which is a constant value 
        that depends on solar wind parameters. All of these functions are independent of beta and gamma. 
        
        Parameters
        ----------
        dipole - dipole tilt angle in radians. 
        pd
        pm
        bz
        
        Returns
        -------
        All coefficients are attached to self. 
        '''

        # Get a coefficients first. 
        a = np.array([12.544, -0.194, 0.305, 0.0573, 2.178, 0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0431, -0.00763, -0.210, 0.0405, -4.430, -0.636, -2.600, 0.832, -5.328, 1.103, -0.907, 1.450])
        #self.a = a

        # Get beta coefficients - renamed delta. 
        beta_c = np.array([a[6] + a[7]*((np.exp(a[8]*bz) - 1)/(np.exp(a[9]*bz) + 1)), a[10], a[11] + a[12]*dipole, a[13]])
         
        # Get cn and cs coefficients (equal). 
        c = a[14]*(pd+pm)**a[15]

        # Get d coefficients. 
        dn = (a[16] + a[17]*dipole + a[18]*dipole**2)
        ds = (a[16] - a[17]*dipole + a[18]*dipole**2)
        
        # Get theta-n and theta-s coefficients.
        theta_n = a[19] + a[20]*dipole
        theta_s = a[19] - a[20]*dipole

        # Get the unscaled subsolar magnetopause radius. 
        r0_lin = 12.544*((pd+pm)**-0.194)*(1 + 0.305*((np.exp(0.0573*bz) -1 )/(np.exp(2.178*bz) + 1)))
        
        return a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin 
        
def get_lin_coeffs_simple(dipole):
        '''This gets the Lin coefficients for the simplified CMEM model that is solar wind independent. 
        Therefore, it does not need to return all of the values. 
        
        Parameters
        ----------
        dipole - dipole tilt angle in radians. 
        
        '''
        
        # Get a coefficients first. 
        a = np.array([12.544, -0.194, 0.305, 0.0573, 2.178, 0.0571, -0.999, 16.473, 0.00152, 0.382, 0.0431, -0.00763, -0.210, 0.0405, -4.430, -0.636, -2.600, 0.832, -5.328, 1.103, -0.907, 1.450])
        
        # Get d coefficients. 
        dn = (a[16] + a[17]*dipole + a[18]*dipole**2)
        ds = (a[16] - a[17]*dipole + a[18]*dipole**2)
        
        # Get theta-n and theta-s coefficients.
        theta_n = a[19] + a[20]*dipole
        theta_s = a[19] - a[20]*dipole
        
        return a, dn, ds, theta_n, theta_s    
        
        
#EMISSIVITY FUNCTIONS.
######################
        
        
def get_model_func(current_model):
        '''This will select the correct function for the desired model. '''
        
        if current_model == "jorg":
            def jorg_func(r, theta, phi, mp, bs, A1, A2, B, alpha, beta, ay_mp, az_mp, ay_bs, az_bs):
               
                '''This is the model from the Jorgensen paper. 
        
                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                mp - subsolar magnetopause distance parameter
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                ay_mp - ay magnetopause flaring parameter
                az_mp - az magnetopause flaring parameter
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''

                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = shue_func(theta, phi, mp, ay_mp, az_mp)
                rbs = shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (A1 + B*((np.sin(theta[r2]))**8))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
        
                return eta
            return jorg_func
             
        elif current_model == "cmem":
            def cmem_func(r, theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs):
                '''
                This is the CMEM model, which will use the lin model to work out 
                the magnetopause location instead of the shue model. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                '''
            
                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3)
                rbs = shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)

                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = A1*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-alpha-(beta*(np.sin(theta[r2]))**2)))
                eta[r3] = A2*((r[r3]/10)**(-3))
                
                return eta
            return cmem_func

        elif current_model == "cmem_linear":
            def cmem_func_linear(r, theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, bs, A1, A2, B, alpha, beta, p1, p2, p3, ay_bs, az_bs, delta=0.5):
                '''
                This is the CMEM model, which will use the lin model to work out 
                the magnetopause location instead of the shue model. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                bs - subsolar bowshock distance parameter
                A1 - parameter
                A2 - parameter
                B - parameter
                alpha - parameter
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                delta - boundary width. def = 0.5 
                
                '''
            
                eta = np.zeros(r.shape)

                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3)
                rbs = shue_func(theta, phi, bs, ay_bs, az_bs)

                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp-delta)
                r2 = np.where((r >= rmp-delta) & (r < rmp))
                r3 = np.where((r >= rmp) & (r < rbs))
                r4 = np.where((r >= rbs) & (r < rbs+delta))
                r5 = np.where(r >= rbs+delta)

                #Calculate eta along the boundaries. 
                eta_rmp = A1*(np.exp(-B*(theta[r2]/2.)**4))*((rmp[r2]/10)**(-alpha-(beta*(np.cos(theta[r2]))**2)))
                
                eta_rbs = A1*(np.exp(-B*(theta[r4]/2.)**4))*((rbs[r4]/10)**(-alpha-(beta*(np.cos(theta[r4]))**2)))
                
                eta_rbs_delta = A2*(((rbs[r4]+delta)/10)**(-3)) 
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (eta_rmp/delta)*(r[r2]-rmp[r2]) + eta_rmp
                eta[r3] = A1*(np.exp(-B*(theta[r3]/2.)**4))*((r[r3]/10)**(-alpha-(beta*(np.sin(theta[r3]))**2)))
                eta[r4] = ((eta_rbs_delta-eta_rbs)/delta)*(r[r4]-rbs[r4]) + eta_rbs 
                eta[r5] = A2*((r[r5]/10)**(-3))
                
                return eta
            return cmem_func_linear
                   
        elif current_model == "cmem2a":
            def cmem2a_func(r, theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, dr, A1, A2, B, dbeta, beta, p1, p2, p3, ay_bs, az_bs):
                '''
                This is the CMEM model, but it is adapted to use constraints. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                dr - distance from MP to BS. (replaced bs in CMEM)
                A1 - parameter
                A2 - parameter
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                
                '''
            
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3)
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0(a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_adapted_func(theta, phi, rmp0, dr, ay_bs, az_bs)
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = A1*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-(dbeta+(beta*(np.cos(theta[r2]))**2))))
                eta[r3] = A2*((r[r3]/10)**(-3))        
            
                return eta
            return cmem2a_func
        
        elif current_model == "cmem2b":
            def cmem2b_func(r, theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, dr, A1, A2, B, dbeta, beta, p1, p2, p3, ay_bs, az_bs, k=50, sigmoid='exp'):
                '''
                This is the CMEM model, which is adapted to use constraints 
                and also sigmoid functions to turn it into one continuous model. #

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin - Lin coefficients in model. 
                p0 - scaling factor on the subsolar magnetopause parameter 
                dr - distance from MP to BS. (replaced bs in CMEM)
                A1 - parameter
                A2 - parameter
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - scaling factor on magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                ay_bs - ay bowshock flaring parameter
                az_bs - az bowshock flaring parameter
                k - boundary width. Assumed equal at MP and BS. Should be provided, not a free parameter at the moment or the fitting will slow down. Default value is 50 if you can't work it out from a cube. Value will be different depending on the sigmoid function chosen. 
                sigmoid - 'exp' or 'arctan' (NOT RECOMMENDED)
                
                
                '''
            
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func(theta, phi, a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p1, p2, p3)
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0(a, beta_c, c, dn, ds, theta_n, theta_s, r0_lin, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_adapted_func(theta, phi, rmp0, dr, ay_bs, az_bs)
        
                #Put new continuous replacement here. 
                #Magnetosheath.
                msheath = A1*(np.exp(-B*(theta/2.)**4))*((r/10)**(-(dbeta+(beta*(np.cos(theta))**2))))
                
                swind = A2*((r/10)**(-3))
                
                if sigmoid == 'exp': 
                    print ('use exp sigmoid')
                    #The warnings filter is included at the top 
                    #for these exponential functions.  
                    s1 = 1/(1+np.exp(-k*(r-rmp))) 
                    s2 = 1/(1+np.exp(k*(r-rbs))) 
                    s3 = 1/(1+np.exp(-k*(r-rbs)))
                    
                elif sigmoid == 'arctan':
                    print ('use arctan sigmoid')
                    print ('DO NOT USE. DOES NOT DECREASE FAST ENOUGH.')
                    print ('YOU WILL GET EMISSION AT LOW VALUES OF R!') 
                    s1 = (np.arctan(k*(r-rmp)) + np.pi/2)/np.pi
                    s2 = (np.arctan(k*(rbs-r)) + np.pi/2)/np.pi
                    s3 = (np.arctan(k*(r-rbs)) + np.pi/2)/np.pi 
                    print (s1.max(), s2.max(), s3.max())
                else:
                    raise ValueError ("Invalid sigmoid function chosen: 'exp' or 'arctan'")     
                    
                eta = msheath*s1*s2 + swind*s3 
        
                return eta
            return cmem2b_func 
        
        elif current_model == 'cmem2c':
            def cmem2c_func(r, theta, phi, a, dn, ds, theta_n, theta_s, p0, dr, A1, A2, B, dbeta, beta, p1, p2, p3, dp1, scaled=False):
                '''
                This is the CMEM model, but it is adapted to use constraints. 

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, dn, ds, theta_n, theta_s - Simplified Lin coefficients in model. 
                p0 - subsolar magnetopause parameter 
                dr - distance from MP to BS. (replaced bs in CMEM)
                A1 - parameter
                A2 - parameter
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                dp1 - Extra flaring on bowshock
                scaled - boolean to use scaled A1 and A2 parameters by 10,000. def = False. 
                
                '''
            
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, p0, p1, p2, p3)
                
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0_simple(a, dn, ds, theta_n, theta_s, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_func_simple(theta, rmp0, dr, p1, dp1)
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                if scaled:
                    eta[r2] = A1*0.00001*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-(dbeta+(beta*(np.cos(theta[r2]))**2))))
                    eta[r3] = A2*0.00001*((r[r3]/10)**(-3)) 
                else:
                    eta[r2] = A1*(np.exp(-B*(theta[r2]/2.)**4))*((r[r2]/10)**(-(dbeta+(beta*(np.cos(theta[r2]))**2))))
                    eta[r3] = A2*((r[r3]/10)**(-3))        
            
                return eta
            return cmem2c_func
        
        elif current_model == 'cmem2d':
            def cmem2d_func(r, theta, phi, a, dn, ds, theta_n, theta_s, p0, dr, A1, A2, B, dbeta, beta, p1, p2, p3, dp1, scaled=False):
                '''
                This is the CMEM 2c model, but it takes absolute 
                values of parameters to avoid using constraints.  

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, dn, ds, theta_n, theta_s - Simplified Lin coefficients in model. 
                p0 - subsolar magnetopause parameter 
                dr - distance from MP to BS. (replaced bs in CMEM)
                A1 - parameter
                A2 - parameter
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                dp1 - Extra flaring on bowshock
                scaled - boolean to use scaled A1 and A2 parameters by 10,000. def = False. 
                
                '''
            
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, abs(p0), abs(p1), abs(p2), abs(p3))
                
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0_simple(a, dn, ds, theta_n, theta_s, abs(p0), abs(p2), abs(p3))
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_func_simple(theta, rmp0, abs(dr), abs(p1), abs(dp1))
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp)
                r2 = np.where((r >= rmp) & (r < rbs))
                r3 = np.where(r >= rbs)
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                if scaled:
                    eta[r2] = abs(A1)*0.00001*(np.exp(-abs(B)*(theta[r2]/2.)**4))*((r[r2]/10)**(-(abs(dbeta)+(abs(beta)*(np.cos(theta[r2]))**2))))
                    eta[r3] = abs(A2)*0.00001*((r[r3]/10)**(-3)) 
                else:
                    eta[r2] = abs(A1)*(np.exp(-abs(B)*(theta[r2]/2.)**4))*((r[r2]/10)**(-(abs(dbeta)+(abs(beta)*(np.cos(theta[r2]))**2))))
                    eta[r3] = abs(A2)*((r[r3]/10)**(-3))        
            
                return eta
            return cmem2d_func

        elif current_model == 'cmem2e':
            def cmem2e_func(r, theta, phi, a, dn, ds, theta_n, theta_s, p0, dr, A1, A2, B, dbeta, beta, p1, p2, p3, dp1, delta=0.5):
                '''
                This is the CMEM 2c model, but with linear sections in the magnetopause and bowshock instead of discontinuous drops.   

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, dn, ds, theta_n, theta_s - Simplified Lin coefficients in model. 
                p0 - subsolar magnetopause parameter 
                dr - distance from MP to BS. (replaced bs in CMEM)
                A1 - parameter
                A2 - parameter
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                dp1 - Extra flaring on bowshock
                delta - boundary width to use. def = 0.5. 
                '''
                
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, p0, p1, p2, p3)
                
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0_simple(a, dn, ds, theta_n, theta_s, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_func_simple(theta, rmp0, dr, p1, dp1)
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp-delta)
                r2 = np.where((r >= rmp-delta) & (r < rmp))
                r3 = np.where((r >= rmp) & (r < rbs))
                r4 = np.where((r >= rbs) & (r < rbs+delta))
                r5 = np.where(r >= rbs+delta)
                
                #Calculate eta along the boundaries. 
                eta_rmp = A1*(np.exp(-B*(theta[r2]/2.)**4))*((rmp[r2]/10)**(-(dbeta+(beta*(np.cos(theta[r2]))**2))))
                
                eta_rbs = A1*(np.exp(-B*(theta[r4]/2.)**4))*((rbs[r4]/10)**(-(dbeta+(beta*(np.cos(theta[r4]))**2))))
                
                eta_rbs_delta = A2*(((rbs[r4]+delta)/10)**(-3)) 
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (eta_rmp/delta)*(r[r2]-rmp[r2]) + eta_rmp
                eta[r3] = A1*(np.exp(-B*(theta[r3]/2.)**4))*((r[r3]/10)**(-(dbeta+(beta*(np.cos(theta[r3]))**2))))
                eta[r4] = ((eta_rbs_delta-eta_rbs)/delta)*(r[r4]-rbs[r4]) + eta_rbs 
                eta[r5] = A2*((r[r5]/10)**(-3))        
            
                return eta
            return cmem2e_func

        elif current_model == 'cmem2f':
            def cmem2f_func(r, theta, phi, a, dn, ds, theta_n, theta_s, p0, A1, B, dbeta, beta, p1, p2, p3, dp1, delta=0.5):
                '''
                This is the CMEM 2c model, but with linear sections in the magnetopause and bowshock instead of discontinuous drops. It also uses the degeneracy between p0 and dr, and A1 and A2, to eliminate dr and A2, instead calculating them from rcmem via p0, and from A1, respectively. These relationships come from CMEM1. See compare_optimised_models.py for the relationships.   

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, dn, ds, theta_n, theta_s - Simplified Lin coefficients in model. 
                p0 - subsolar magnetopause parameter. dr = 0.626rcmem-3.166
                A1 - parameter. Old A2 is calculated from A1. A2 = 0.41A1 
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter
                p1 - magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                dp1 - Extra flaring on bowshock
                delta - boundary width to use. def = 0.5. 
                '''
                
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, p0, p1, p2, p3)
                
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0_simple(a, dn, ds, theta_n, theta_s, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_func_simple_cmem2f(theta, rmp0, p1, dp1)
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp-delta)
                r2 = np.where((r >= rmp-delta) & (r < rmp))
                r3 = np.where((r >= rmp) & (r < rbs))
                r4 = np.where((r >= rbs) & (r < rbs+delta))
                r5 = np.where(r >= rbs+delta)
                
                #Calculate eta along the boundaries. 
                eta_rmp = A1*(np.exp(-B*(theta[r2]/2.)**4))*((rmp[r2]/10)**(-(dbeta+(beta*(np.cos(theta[r2]))**2))))
                
                eta_rbs = A1*(np.exp(-B*(theta[r4]/2.)**4))*((rbs[r4]/10)**(-(dbeta+(beta*(np.cos(theta[r4]))**2))))
                
                eta_rbs_delta = A1*0.41*(((rbs[r4]+delta)/10)**(-3)) 
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (eta_rmp/delta)*(r[r2]-rmp[r2]) + eta_rmp
                eta[r3] = A1*(np.exp(-B*(theta[r3]/2.)**4))*((r[r3]/10)**(-(dbeta+(beta*(np.cos(theta[r3]))**2))))
                eta[r4] = ((eta_rbs_delta-eta_rbs)/delta)*(r[r4]-rbs[r4]) + eta_rbs 
                eta[r5] = A1*0.41*((r[r5]/10)**(-3))      
                #Used relationship with A1 to eliminate A2  
            
                return eta
            return cmem2f_func
 
        elif current_model == 'cmem2g':
            def cmem2g_func(r, theta, phi, a, dn, ds, theta_n, theta_s, p0, A1, B, dbeta, p1, p2, p3, delta=0.5):
                '''
                This is the CMEM 2c model, but with linear sections in the magnetopause and bowshock instead of discontinuous drops. It also uses the degeneracy between p0 and dr, and A1 and A2, to eliminate dr and A2, instead calculating them from rcmem via p0, and from A1, respectively. These relationships come from CMEM1. See compare_optimised_models.py for the relationships. Also fixes beta and dp1.   

                Parameters
                ----------
                r - 3D array of r values.
                theta - 3D array of theta values. 
                phi - 3D array of phi values. 
                a, dn, ds, theta_n, theta_s - Simplified Lin coefficients in model. 
                p0 - subsolar magnetopause parameter. dr = 0.626rcmem-3.166
                A1 - parameter. Old A2 is calculated from A1. A2 = 0.41A1 
                B - parameter
                dbeta - parameter (replaced alpha in CMEM) 
                beta - parameter. Fixed at 2. 
                p1 - magnetopause flaring parameter
                p2 - scaling parameter on magnetopause indentation parameter 
                p3 - scaling parameter on magnetopause indentation parameter 
                dp1 - Extra flaring on bowshock. Fixed at 0.05. 
                delta - boundary width to use. def = 0.5. 
                '''
                
                eta = np.zeros(r.shape)
                
                # Calculate the radii to the magnetopause and bowshock for all 
                # combinations of theta and phi. 
                rmp = lin_scaled_func_simple(theta, phi, a, dn, ds, theta_n, theta_s, p0, p1, p2, p3)
                
                #Get subsolar magnetopause value. 
                rmp0 = get_rmp0_simple(a, dn, ds, theta_n, theta_s, p0, p2, p3)
                
                #Use the adapted function to get the bowshock position. 
                rbs = shue_func_simple_cmem2f(theta, rmp0, p1, 0.05)
        
                # Get indices inside MP, between MP and BS, and outside BS. 
                r1 = np.where(r < rmp-delta)
                r2 = np.where((r >= rmp-delta) & (r < rmp))
                r3 = np.where((r >= rmp) & (r < rbs))
                r4 = np.where((r >= rbs) & (r < rbs+delta))
                r5 = np.where(r >= rbs+delta)
                
                #Calculate eta along the boundaries. 
                eta_rmp = A1*(np.exp(-B*(theta[r2]/2.)**4))*((rmp[r2]/10)**(-(dbeta+(2*(np.cos(theta[r2]))**2))))
                
                eta_rbs = A1*(np.exp(-B*(theta[r4]/2.)**4))*((rbs[r4]/10)**(-(dbeta+(2*(np.cos(theta[r4]))**2))))
                
                eta_rbs_delta = A1*0.41*(((rbs[r4]+delta)/10)**(-3)) 
                
                # Now calculate eta in each region. 
                eta[r1] = 0.0
                eta[r2] = (eta_rmp/delta)*(r[r2]-rmp[r2]) + eta_rmp
                eta[r3] = A1*(np.exp(-B*(theta[r3]/2.)**4))*((r[r3]/10)**(-(dbeta+(2*(np.cos(theta[r3]))**2))))
                eta[r4] = ((eta_rbs_delta-eta_rbs)/delta)*(r[r4]-rbs[r4]) + eta_rbs 
                eta[r5] = A1*0.41*((r[r5]/10)**(-3))      
                #Used relationship with A1 to eliminate A2  
            
                return eta
            return cmem2g_func
                                             
        else:
            raise ValueError("{} not a valid model. 'jorg', 'cmem', 'cmem2a', cmem2b', cmem2c', 'cmem2d', 'cmem2e', 'cmem2f' or 'cmem2g' only atm.".format(current_model))

