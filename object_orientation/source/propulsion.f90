module propulsion_m 
    use koch_m 
    use atmosphere_m 
    use controller_m
    use jsonx_m 
    implicit none  

    contains 

    subroutine propulsion_init(t,j_propulsion) 
        implicit none 
        type(propulsion_t), intent(inout) :: t 
        type(json_value), pointer :: j_propulsion 
        integer :: i 
        character(len=:), allocatable :: temp 

        t%name = j_propulsion%name 
        write(*,*) '    -Initializing Propulsion : ', trim(t%name) 

        call jsonx_get(j_propulsion, 'location[ft]', t%location, 0.0, 3) 
        call jsonx_get(j_propulsion, 'orientation[deg]', t%orientation_eul, 0.0, 3) 
        call jsonx_get(j_propulsion, 'type', t%type) 
        t%orientation_eul = t%orientation_eul * pi / 180.0 
        t%orientation_quat = euler_to_quat(t%orientation_eul) 

        select case(t%type) 
            case("T=f(V)") 
                call jsonx_get(j_propulsion, 'T_coefficients[lbf]', t%T_coeffs, 0.0) 
                call jsonx_get(j_propulsion, 'Ta', t%Ta) 
                t%rotation_delta = 1 
                t%Ixx = 0.0 
            case("propeller_polynomial") 
                call jsonx_get(j_propulsion, 'diameter[ft]', t%diameter) 
                call jsonx_get(j_propulsion, 'Ixx[slug-ft^2]', t%Ixx) 
                call jsonx_get(j_propulsion, 'rotation', temp) 
                t%rotation_delta = 1
                if(trim(temp) == 'LH') t%rotation_delta = -1 
                call jsonx_get(j_propulsion, 'CT(J)', t%CT_J, 0.0) 
                call jsonx_get(j_propulsion, 'CPb(J)', t%CP_J, 0.0) 
                call jsonx_get(j_propulsion, 'CN,alpha(J)', t%CNa_J, 0.0) 
                call jsonx_get(j_propulsion, 'Cn,alpha(J)', t%Cnna_J, 0.0) 
        end select 
    end subroutine propulsion_init 

    function propulsion_get_FMh(t, states, tau) result(ans) 
        implicit none 
        type(propulsion_t), intent(inout) :: t 
        real, intent(in) :: states(21), tau 
        real :: ans(9) 
        real :: Vc(3), Vc_mag, uc(3), vN(3), vN_mag, uN(3) 
        real :: Fc(3), Mc(3) 
        real :: thrust, normal, torque, yaw, hxx, alpha_c 
        real :: Z_dum, T_dum, P_dum, rho, rho0, a_dum, mu_dum, dyp 
        real :: Hz, omega, J

        Vc = quat_base_to_dependent(states(1:3) + cross_product(states(4:6), t%location), t%orientation_quat)
        Vc_mag = sqrt(Vc(1)**2 + Vc(2)**2 + Vc(3)**2) 
        uc = Vc/Vc_mag 
        if(Vc_mag < tol) uc = [1.0, 0.0, 0.0] 
        alpha_c = acos(uc(1))

        vN = -[0.0, uc(2), uc(3)] 
        vN_mag = sqrt(vN(1)**2 + vN(2)**2 + vN(3)**2) 
        uN = vN/ vN_mag 
        if(vN_mag < tol) uN = [0.0, 0.0, 1.0] 

        call std_atm_English(0.0, z_dum, t_dum, p_dum, rho0, a_dum, mu_dum) 
        call std_atm_English(-states(9), z_dum, t_dum, p_dum, rho, a_dum, mu_dum)

        select case(t%type) 
            case("T=f(V)") 
                thrust = tau*calc_polynomial(t%T_coeffs, Vc_mag) * (rho/rho0)**t%Ta 
                normal = 0.0 
                torque = 0.0 
                yaw = 0.0
                hxx = 0.0 
            case("propeller_polynomial") 
                select case(t%units)
                    case("RPM")
                        Hz = tau/60.0 
                        omega = Hz*2*pi 
                    case ("thrust[lbf]")
                        omega = -Vc_mag * t%CT_J(2) + sqrt(Vc_mag**2 * t%CT_J(2)**2 - 4.0 * t%CT_J(1)*(Vc_mag**2*t%CT_J(3) - tau / rho/ t%diameter**2) )
                        omega = omega * pi /t%diameter / t%CT_J(1) 
                        Hz = omega / 2.0 / pi 
                    case("torque[ft-lbf]") 
                        omega = -Vc_mag * t%CP_J(2) + sqrt(Vc_mag**2*t%CP_J(2)**2 - 4.0 * t%CP_J(1)*(Vc_mag**2*t%CP_J(3) - tau *2.0 * pi / rho / t%diameter**3) ) 
                        omega = omega * pi / t%diameter / t%CP_J(1) 
                        Hz = omega / 2.0 / pi 
                end select

                    ! write(*,*) 'omega = ', omega 
                    ! write(*,*) 'J = ', J                 
                    ! write(*,*) 't%CT_J = ', t%CT_J
                    ! write(*,*) 't%CP_J = ', t%CP_J
                    ! write(*,*) 't%CNa_J = ', t%CNa_J
                    ! write(*,*) 't%Cnna_J = ', t%Cnna_J

                    J = 2.0 * pi * Vc_mag/omega/t%diameter 
                    thrust = calc_polynomial(t%CT_J,J)    * rho*(Hz**2)*(t%diameter**4) 
                    torque = calc_polynomial(t%CP_J,J)    * rho*(Hz**3)*(T%diameter**5) / omega 
                    normal = calc_polynomial(t%CNa_J,J)   * rho*(Hz**2)*(t%diameter**4) * alpha_c 
                    yaw    = calc_polynomial(t%Cnna_J,J)  * rho*(Hz**2)*(T%diameter**5) * alpha_c 
                    hxx = t%rotation_delta * t%Ixx * omega 
                    ! write(*,*) 
                    ! write(*,*) 'thrust =', thrust
                    ! write(*,*) 'torque =', torque
                    ! write(*,*) 'normal =', normal
                    ! write(*,*) 'yaw    =', yaw   
                    ! write(*,*) 'alpha_c = ', alpha_c 
                    ! write(*,*) 'hxx    =', hxx
            end select 
        ! write(*,*) 'uN = ', uN

        Fc = [thrust, 0.0, 0.0] + Normal*uN 
        Mc = -real(t%rotation_delta)*([torque, 0.0, 0.0] + yaw*uN) 

        ans(1:3) = quat_dependent_to_base(Fc, t%orientation_quat) 
        ans(4:6) = quat_dependent_to_base(Mc, t%orientation_quat) + cross_product(t%location, ans(1:3)) 
        ans(7:9) = quat_dependent_to_base([hxx, 0.0, 0.0], t%orientation_quat)
        ! write(*,*)
        ! write(*,*) t%name 
        ! write(*,*) 'Fc,','Mc',',',Fc(1),',',Fc(2),',',Fc(3),',',Mc(1),',',Mc(2),',',Mc(3)        
        ! write(*,*) 'Fb',',','Mb',',',ans(1),',',ans(2),',',ans(3),',',ans(4),',',ans(5),',',ans(6)
        ! write(*,*) 'hc',',','hb',',',hxx,',',0.0,',',0.0,',',ans(7),',',ans(8),',',ans(9)
    end function propulsion_get_FMh

    function calc_polynomial(coeffs, var) result(ans) 
        implicit none 
        real, intent(in) :: coeffs(:), var
        real :: ans 
        integer :: i 
        ans = 0.0 
        do i=1,size(coeffs) 
            ans = ans + coeffs(i)*var**(i-1) 
        end do 
    end function calc_polynomial 

    function electric_propulsion(tau, Vc, gamma, tolerance) result(Nr)
        real, intent(in) :: tau, Vc, gamma, tolerance 
        real :: Nr, error 
        real :: lw, ls, lr, lm, Im, eta_c, B, C, Em, Nm, Ns, error 
        real :: Eb, Ib 
        real :: Gm, eta_g, Kv, CI, Im0, Rc, Eb0, Rb, Rm, Rc

        Nr = 1000 
        error = 1.0 
        do while (error > tolerance)
            ls = lr 
            lm = ls/(eta_g * Gm) 
            Im = lm * Kv/CI + Im0 
            eta_c = 1 - 0.078*(1-tau) 
            B = 2*Im*Rc - tau*Eb0 + tau**2/eta_c *Im * Rb 
            C = Im**2 * Rc**2 - tau*Eb0*Im*Rc 
            Em = 0.5*(-B + sqrt(B**2 -4C))
            Nm = Kv * (Em - Im*Rm)
            Ns = Nm/Gm 
            error = Ns - Nr 
            Nr = Nr + gamma*(Ns-Nr)
        end do 
        Eb = (Nm/Kv + Im*(Rm+Rc))/tau 
        Ib = (Eb0 - Eb)/Rb 


end module propulsion_m



            
