/*
 *
 * Copyright (c) 2019 Ewoud Smeur and Evghenii Volodscoi
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * This control algorithm is the Incremental Nonlinear Dynamic Inversion (INDI)
 * controller to control the position of the Crazyflie. It can be seen as an extension
 * (outer loop) to the already existing INDI attitude controller (inner loop).
 * 
 * The control algorithm was implemented according to the publication in the
 * journal od Control Engineering Practice: Cascaded Incremental Nonlinear Dynamic 
 * Inversion for MAV Disturbance Rejection
 * https://doi.org/10.1016/j.conengprac.2018.01.003
 */


#include "position_controller_indi.h"
#include "math3d.h"
#include "pm.h"
#include "range.h"
#include "power_distribution.h"
#include "motors.h"

// Position controller gains
float K_xi_x = 1.0f;
float K_xi_y = 1.0f;
float K_xi_z = 1.0f;
// Velocity controller gains
float K_dxi_x = 5.0f;
float K_dxi_y = 5.0f;
float K_dxi_z = 5.0f;
// Thrust mapping parameter
// float K_thr = 0.00024730f;
float K_thr = 0.0003;


float mass = 0.033;
// float mass = 0.037;

float Thrust_0 = 0.0f;
float Thrust_0_N = 0.0f;
float motor_norm[] = {0, 0, 0, 0};
float motor_norm_test[] = {0,0,0,0};

float g11; float g12; float g13; float g21; float g22; float g23; float g31; float g32; float g33; 
float g14; float g24; float g34;

static float posS_x, posS_y, posS_z;			// Current position
static float velS_x, velS_y, velS_z;			// Current velocity

// Reference values
static struct Vectr positionRef; 
static struct Vectr velocityRef;

// Ceiling effect Estimation
float ceiling_dist = 1.0f;
float T_Ratio_RBF  = 1.0f; //variable to store ranger measurement 
float RBF_der_value = 1.0f;
float T_Ratio_meas;
float K_ctr_eff = 0.002;


/** Stucture for the RBF function */
static struct RBF_form_s RBF_t = {.centers_s = {0, 0.005, 0.015, 0.03, 0.06, 0.1}
								, .kernel_width_s = {15000.0, 10000.0, 5000.0, 2000.0, 1250.0, 1000.0}
								// , .weights_s = {1.43142386f,  0.15021046f, -0.09754395f ,-0.07601899f ,-0.26694495f, -0.30407099f, -0.37406497f}};
								// , .weights_s = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}};
								, .weights_s = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}};


/** Structure for Ordinary least square regression (OLS) */
static struct OLS_s OLS_t = {.reg_rows = 21
							, .reg_cols = 7
							, .reg_matrix_s = {0.0f}
						 	, .data_points_x_s = {0.0f, 0.005f, 0.01f, 0.015f, 0.02f, 0.025f,  0.030f, 0.035f,  0.040f, 0.045f, 0.050f, 0.055f, 0.060f, 0.065f, 0.070f, 0.075f,  0.080f, 0.085f, 0.090f, 0.095f, 0.10f}
							, .data_points_y_s = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}
							};


/** calculate one RBF function value, needed in regression matrix, thus calculated withouth the weight factor*/
float RBF_element(struct RBF_form_s *RBF, float x, int element){
	float element_result;
	element_result = expf(-RBF->kernel_width_s[element]*powf((x - RBF->centers_s[element]), 2));	
	return element_result;
};

/** calculated the value at given point of the RBF function */
float RBF_result(struct RBF_form_s *RBF, float x) {
	float temp;
	float total = 0;
	//Calculate value of the 6 RBF hills and add together
	for (int i = 0; i < 6; i++){
		temp = RBF->weights_s[i+1]*expf(-RBF->kernel_width_s[i]*powf((x - RBF->centers_s[i]), 2));	
		total += temp;
	}

	//Add the constant times its weight
	total += RBF->weights_s[0]*1;
	
	return total;
};

/** calculated the value at given point of the derivative of the RBF function */
float RBF_der_result(struct RBF_form_s *RBF, float x) {
	float temp = 0;
	
	for (int i = 0; i < 6; i++){
		temp += -2*RBF->kernel_width_s[i]*(x-RBF->centers_s[i])*RBF->weights_s[i+1]*expf(-RBF->kernel_width_s[i]*powf((x - RBF->centers_s[i]), 2));
	}

	return temp;
};

/** Calculate the pseudo inverse of the refression matrix*/
arm_matrix_instance_f32 Pseudo_Inverse(struct OLS_s *OLS){
	arm_matrix_instance_f32 A; //regression matrix instance
	arm_mat_init_f32(&A, OLS->reg_rows, OLS->reg_cols, (float32_t *)OLS->reg_matrix_s);

	float32_t AT_buf[7*21];
	arm_matrix_instance_f32 AT; //transpose regression matrix instance
	arm_mat_init_f32(&AT, 7, 21, (float32_t *)AT_buf);
	arm_mat_trans_f32(&A, &AT);

	float32_t AT_A_buf[7*7];
	arm_matrix_instance_f32 AT_A; //Multiplication of AT and A 
	arm_mat_init_f32(&AT_A, 7, 7, (float32_t *)AT_A_buf);
	arm_mat_mult_f32(&AT, &A, &AT_A);

	float32_t AT_A_inv_buf[7*7];
	arm_matrix_instance_f32 AT_A_inv; //Inverse of AT times A
	arm_mat_init_f32(&AT_A_inv, 7, 7, (float32_t *)AT_A_inv_buf);
	arm_mat_inverse_f32(&AT_A, &AT_A_inv);

	float32_t pseudo_inv_buf[7*21];
	arm_matrix_instance_f32 pseudo_inv;
	arm_mat_init_f32(&pseudo_inv, 7, 21, (float32_t *)pseudo_inv_buf);
	arm_mat_mult_f32(&AT_A_inv, &AT, &pseudo_inv); 

	return pseudo_inv;
};

/** Calculated the weights with OLS*/
void OLS_weights(struct RBF_form_s *RBF, struct OLS_s *OLS){
	arm_matrix_instance_f32 RBF_data_value; 
	arm_mat_init_f32(&RBF_data_value, 21, 1, (float32_t *)OLS->data_points_y_s);

	float32_t new_weights_buf[7]; 
	arm_matrix_instance_f32 new_weights; 
	arm_mat_init_f32(&new_weights, 7, 1, (float32_t*)new_weights_buf);

	arm_matrix_instance_f32 pseudo_inv;
	pseudo_inv = Pseudo_Inverse(OLS);

	arm_mat_mult_f32(&pseudo_inv, &RBF_data_value, &new_weights);
	
	for (int i = 0; i < 7; i++){
		RBF->weights_s[i] = new_weights.pData[i];
	}
}

/** Create the regression matrix, only if new data point is accepted*/
void Create_reg_mat(struct OLS_s *OLS, struct RBF_form_s *RBF){
	int j = 0;
	int k = 0;

	for (int i = 0; i < 7*21; i++){

		if (i != 0 && i%7 == 0){
			k++;	
		}

		if (i == 0 || i%7 == 0){
			OLS->reg_matrix_s[i] = 1.0f;
			j = 0;
		}

		else{
			OLS->reg_matrix_s[i] = RBF_element(RBF,  OLS->data_points_x_s[k], j);
			j++;
		};
	};
};

void LMS_weight(struct RBF_form_s *RBF, float RBF_value, float RBF_meas_value, float data_point){
	for (int i = 1; i < 7; i++){
		RBF->weights_s[i] = RBF->weights_s[i] + 0.05f*(RBF_meas_value - RBF_value)*RBF_element(RBF, data_point, i);
	}

	RBF->weights_s[0] = RBF->weights_s[0] + 0.05f*(RBF_meas_value - RBF_value);
}



static struct IndiOuterVariables indiOuter = {
	.filt_cutoff = POSITION_INDI_FILT_CUTOFF,
	.act_dyn_posINDI = STABILIZATION_INDI_ACT_DYN_P
};


void position_indi_init_filters(void)
{
	// tau = 1/(2*pi*Fc)
	float tau = 1.0f / (2.0f * M_PI_F * indiOuter.filt_cutoff);
	float tau_axis[3] = {tau, tau, tau};
	float sample_time = 1.0f / ATTITUDE_RATE;
	// Filtering of linear acceleration, attitude and thrust 
	for (int8_t i = 0; i < 3; i++) {
		init_butterworth_2_low_pass(&indiOuter.ddxi[i], tau_axis[i], sample_time, 0.0f);
		init_butterworth_2_low_pass(&indiOuter.ang[i], tau_axis[i], sample_time, 0.0f);
		init_butterworth_2_low_pass(&indiOuter.thr[i], tau_axis[i], sample_time, 0.0f);
	}
}

// Linear acceleration filter
static inline void filter_ddxi(Butterworth2LowPass *filter, struct Vectr *old_values, struct Vectr *new_values)
{
	new_values->x = update_butterworth_2_low_pass(&filter[0], old_values->x);
	new_values->y = update_butterworth_2_low_pass(&filter[1], old_values->y);
	new_values->z = update_butterworth_2_low_pass(&filter[2], old_values->z);
}

// Attitude filter
static inline void filter_ang(Butterworth2LowPass *filter, struct Angles *old_values, struct Angles *new_values)
{
	new_values->phi = update_butterworth_2_low_pass(&filter[0], old_values->phi);
	new_values->theta = update_butterworth_2_low_pass(&filter[1], old_values->theta);
	new_values->psi = update_butterworth_2_low_pass(&filter[2], old_values->psi);
}

// Thrust filter
static inline void filter_thrust(Butterworth2LowPass *filter, float *old_thrust, float *new_thrust) 
{
	*new_thrust = update_butterworth_2_low_pass(&filter[0], *old_thrust);
}


// Computes transformation matrix from body frame (index B) into NED frame (index O)
void m_ob(struct Angles att, float matrix[3][3]) {

	matrix[0][0] = cosf(att.theta)*cosf(att.psi);
	matrix[0][1] = sinf(att.phi)*sinf(att.theta)*cosf(att.psi) - cosf(att.phi)*sinf(att.psi); 
	matrix[0][2] = cosf(att.phi)*sinf(att.theta)*cosf(att.psi) + sinf(att.phi)*sinf(att.psi);
	matrix[1][0] = cosf(att.theta)*sinf(att.psi);
	matrix[1][1] = sinf(att.phi)*sinf(att.theta)*sinf(att.psi) + cosf(att.phi)*cosf(att.psi);
	matrix[1][2] = cosf(att.phi)*sinf(att.theta)*sinf(att.psi) - sinf(att.phi)*cosf(att.psi);
	matrix[2][0] = -sinf(att.theta);
	matrix[2][1] = sinf(att.phi)*cosf(att.theta);
	matrix[2][2] = cosf(att.phi)*cosf(att.theta);
}


void positionControllerINDIInit(void)
{
	// Re-initialize filters
	position_indi_init_filters();
}


void positionControllerINDI(const sensorData_t *sensors,
                            setpoint_t *setpoint,
                            const state_t *state, 
                            vector_t *refOuterINDI){ 

	// Read states (position, velocity)
	posS_x = state->position.x;
	posS_y = -state->position.y;
	posS_z = -state->position.z;
	velS_x = state->velocity.x;
	velS_y = -state->velocity.y;
	velS_z = -state->velocity.z;

	// Read in velocity setpoints
    velocityRef.x = setpoint->velocity.x;
	velocityRef.y = -setpoint->velocity.y;
    velocityRef.z = -setpoint->velocity.z;

	// Position controller (Proportional)
	if (setpoint->mode.x == modeAbs) {
		positionRef.x = setpoint->position.x;
		velocityRef.x = K_xi_x*(positionRef.x - posS_x);
	}
	if (setpoint->mode.y == modeAbs) {
		positionRef.y = -setpoint->position.y;
		velocityRef.y = K_xi_y*(positionRef.y - posS_y);
	}
	if (setpoint->mode.z == modeAbs) {
		positionRef.z = -setpoint->position.z;
		velocityRef.z = K_xi_z*(positionRef.z - posS_z);
	}

	// Velocity controller (Proportional)
	indiOuter.linear_accel_ref.x = K_dxi_x*(velocityRef.x - velS_x);
	indiOuter.linear_accel_ref.y = K_dxi_y*(velocityRef.y - velS_y);
	indiOuter.linear_accel_ref.z = K_dxi_z*(velocityRef.z - velS_z); 

	// Acceleration controller (INDI)
	// Read lin. acceleration (Body-fixed) obtained from sensors CHECKED
	indiOuter.linear_accel_s.x = (sensors->acc.x)*9.81f;
	indiOuter.linear_accel_s.y = (-sensors->acc.y)*9.81f;
	indiOuter.linear_accel_s.z = (-sensors->acc.z)*9.81f;

	// Filter lin. acceleration 
	filter_ddxi(indiOuter.ddxi, &indiOuter.linear_accel_s, &indiOuter.linear_accel_f);

	// Obtain actual attitude values (in rad)
	indiOuter.attitude_s.phi = radians(state->attitude.roll); 
	indiOuter.attitude_s.theta = radians(state->attitude.pitch);
	indiOuter.attitude_s.psi = -radians(state->attitude.yaw);
	filter_ang(indiOuter.ang, &indiOuter.attitude_s, &indiOuter.attitude_f);


	// Actual attitude (in rad)
	struct Angles att = {
		.phi = indiOuter.attitude_f.phi,
		.theta = indiOuter.attitude_f.theta,
		.psi = indiOuter.attitude_f.psi,
	};

	// Compute transformation matrix from body frame (index B) into NED frame (index O)
	float M_OB[3][3] = {0};
	m_ob(att, M_OB);

	// Transform lin. acceleration in NED (add gravity to the z-component)
	indiOuter.linear_accel_ft.x = M_OB[0][0]*indiOuter.linear_accel_f.x + M_OB[0][1]*indiOuter.linear_accel_f.y + M_OB[0][2]*indiOuter.linear_accel_f.z;
	indiOuter.linear_accel_ft.y = M_OB[1][0]*indiOuter.linear_accel_f.x + M_OB[1][1]*indiOuter.linear_accel_f.y + M_OB[1][2]*indiOuter.linear_accel_f.z;
	indiOuter.linear_accel_ft.z = M_OB[2][0]*indiOuter.linear_accel_f.x + M_OB[2][1]*indiOuter.linear_accel_f.y + M_OB[2][2]*indiOuter.linear_accel_f.z + 9.81f; 

	// Compute lin. acceleration error
	indiOuter.linear_accel_err.x = indiOuter.linear_accel_ref.x - indiOuter.linear_accel_ft.x;
	indiOuter.linear_accel_err.y = indiOuter.linear_accel_ref.y - indiOuter.linear_accel_ft.y;
	indiOuter.linear_accel_err.z = indiOuter.linear_accel_ref.z - indiOuter.linear_accel_ft.z;

	// Elements of the G matrix (see publication for more information) 
	// ("-" because T points in neg. z-direction, "*9.81" because T/m=a=g, 
	// negative psi to account for wrong coordinate frame in the implementation of the inner loop)
	float Battery_norm = pmGetBatteryVoltage()/NORM_BATTERY;

	motor_norm_test[0] = (motorsCompensateBatteryVoltage(motorPower.m1)/NORM_THRUST);
	motor_norm_test[1] = (motorsCompensateBatteryVoltage(motorPower.m2)/NORM_THRUST);
	motor_norm_test[2] = (motorsCompensateBatteryVoltage(motorPower.m3)/NORM_THRUST);
	motor_norm_test[3] = (motorsCompensateBatteryVoltage(motorPower.m4)/NORM_THRUST);

	motor_norm[0] = 11.093358483549203f - 39.08104165843915f * (motorPower.m1/NORM_THRUST) - 9.525647087583181f * Battery_norm + 20.573302305476638f * (motorPower.m1/NORM_THRUST) * (motorPower.m1/NORM_THRUST) + 38.42885066644033f * (motorPower.m1/NORM_THRUST) * Battery_norm;
	motor_norm[1] = 11.093358483549203f - 39.08104165843915f * (motorPower.m2/NORM_THRUST) - 9.525647087583181f * Battery_norm + 20.573302305476638f * (motorPower.m2/NORM_THRUST) * (motorPower.m2/NORM_THRUST) + 38.42885066644033f * (motorPower.m2/NORM_THRUST) * Battery_norm;
	motor_norm[2] = 11.093358483549203f - 39.08104165843915f * (motorPower.m3/NORM_THRUST) - 9.525647087583181f * Battery_norm + 20.573302305476638f * (motorPower.m3/NORM_THRUST) * (motorPower.m3/NORM_THRUST) + 38.42885066644033f * (motorPower.m3/NORM_THRUST) * Battery_norm;
	motor_norm[3] = 11.093358483549203f - 39.08104165843915f * (motorPower.m4/NORM_THRUST) - 9.525647087583181f * Battery_norm + 20.573302305476638f * (motorPower.m4/NORM_THRUST) * (motorPower.m4/NORM_THRUST) + 38.42885066644033f * (motorPower.m4/NORM_THRUST) * Battery_norm;

	float sum;
	int loop;
   	sum = 0;
   
    for(loop = 3; loop >= 0; loop--) {
		sum = sum + motor_norm[loop];      
  	}

	Thrust_0_N = Thrust_0_N + indiOuter.act_dyn_posINDI*(sum - Thrust_0_N);
	// Thrust_0_N = sum;

	Thrust_0 = ((Thrust_0_N/1000.0f)*-9.81f)/mass; 

	//Thrust_0 cannot become zero due to matrix inversion.
	if (Thrust_0 > -(9.81f/2.0f)){
		Thrust_0 = -9.81f/2.0f;
	};

	//RBF calculation
	//Step 1 get distance and Thrust ratio measurement
	logVarId_t idUp = logGetVarId("range", "up");
	uint16_t left = logGetUint(idUp);
	ceiling_dist = left/1000.0f;
	// ceiling_dist = 0.01;

	bool T_ratio_accept;

	if (indiOuter.linear_accel_s.z < 0){
		T_Ratio_meas = fabs(indiOuter.linear_accel_s.z / Thrust_0); //both are in m/s^2 (for ratio units do not matter) KAN 0 Zijn 
		T_ratio_accept = true; 
	}
	else{
		T_ratio_accept = false;
	}

	// ceiling_dist = 0.5;
	// T_Ratio_meas = 1.2;
	// T_ratio_accept = true;

	// Step 2 see if the measurement can be accepted as new data point
	// bool data_accept = false;

 	if (ceiling_dist <= 0.1f && T_Ratio_meas < 2.5f && T_ratio_accept == true){ //check if close enough at ceiling, and cap the max tratio to count for weird noise data points
		
		float data_dist_smallest = fabs(OLS_t.data_points_x_s[0] - ceiling_dist);
		int idx_data_dist_smallest = 0;

		for (int i = 1; i < 21; i++){
			float data_dist = fabs(OLS_t.data_points_x_s[i] - ceiling_dist);
			
			if (data_dist < data_dist_smallest){
				data_dist_smallest = data_dist;
				idx_data_dist_smallest = i;
			}
		}

		if (data_dist_smallest < 0.0025f && (T_Ratio_meas >  OLS_t.data_points_y_s[idx_data_dist_smallest]*1.1f || T_Ratio_meas < OLS_t.data_points_y_s[idx_data_dist_smallest]*0.9f)){
			// data_accept = true; 

			for (int i = 0; i < idx_data_dist_smallest + 1; i++){
				OLS_t.data_points_y_s[i] = T_Ratio_meas;
			}

			LMS_weight(&RBF_t, RBF_result(&RBF_t, OLS_t.data_points_x_s[idx_data_dist_smallest]), T_Ratio_meas, OLS_t.data_points_x_s[idx_data_dist_smallest]);
		}
	}

	// // //Step 3 update the regression matrix with new data point
	// if (data_accept){
	// // 	Create_reg_mat(&OLS_t, &RBF_t);
	
	// //Step 4 Calculate new weights with pseudoinverse
	// 	// OLS_weights(&RBF_t, &OLS_t);
	// }

	//Step 5 Calculate the thrust ratio from the rbf and use that in the control effectiveness matrix
	if (ceiling_dist <= 0.1f){
		T_Ratio_RBF = RBF_result(&RBF_t, ceiling_dist); // the RBF ratio is in the control effectiveness when deriving to phi,theta,T
		RBF_der_value = RBF_der_result(&RBF_t, ceiling_dist); //the derivative of the RBF is in the control effectiveness when deriving to z

		// if (T_Ratio_RBF > 3){
		// 	T_Ratio_RBF = 1.5;
		// }

		if (T_Ratio_RBF < 1){
			T_Ratio_RBF = 1;
		}

	}

	//END RBF CALC

	if (ceiling_dist <= 0.1f){
		g11 = (cosf(att.phi)*sinf(att.psi) - sinf(att.phi)*sinf(att.theta)*cosf(att.psi))*Thrust_0*T_Ratio_RBF; //(-9.81f); 
		g12 = (cosf(att.phi)*cosf(att.theta)*cosf(att.psi))*Thrust_0*T_Ratio_RBF; //(-9.81f); 						
		g13 = (sinf(att.phi)*sinf(att.psi) + cosf(att.phi)*sinf(att.theta)*cosf(att.psi))*T_Ratio_RBF;
		g21 = (-cosf(att.phi)*cosf(att.psi) - sinf(att.phi)*sinf(att.theta)*sinf(att.psi))*Thrust_0*T_Ratio_RBF; //(-9.81f);
		g22 = (cosf(att.phi)*cosf(att.theta)*sinf(att.psi))*Thrust_0*T_Ratio_RBF; //(-9.81f);
		g23 = (-sinf(att.phi)*cosf(att.psi) + cosf(att.phi)*sinf(att.theta)*sinf(att.psi))*T_Ratio_RBF;
		g31 = (-sinf(att.phi)*cosf(att.theta))*Thrust_0*T_Ratio_RBF; //(-9.81f);
		g32 = (-cosf(att.phi)*sinf(att.theta))*Thrust_0*T_Ratio_RBF; //(-9.81f);
		g33 = (cosf(att.phi)*cosf(att.theta))*T_Ratio_RBF;

		g14 = (sinf(att.phi)*sinf(att.psi) + cosf(att.phi)*sinf(att.theta)*cosf(att.psi))*Thrust_0*RBF_der_value;
		g24 = (-sinf(att.phi)*cosf(att.psi) + cosf(att.phi)*sinf(att.theta)*sinf(att.psi))*Thrust_0*RBF_der_value;
		g34 = (cosf(att.phi)*cosf(att.theta))*Thrust_0*RBF_der_value;
	}
	else {
		g11 = (cosf(att.phi)*sinf(att.psi) - sinf(att.phi)*sinf(att.theta)*cosf(att.psi))*Thrust_0; //(-9.81f); 
		g12 = (cosf(att.phi)*cosf(att.theta)*cosf(att.psi))*Thrust_0; //(-9.81f); 						
		g13 = (sinf(att.phi)*sinf(att.psi) + cosf(att.phi)*sinf(att.theta)*cosf(att.psi));
		g21 = (-cosf(att.phi)*cosf(att.psi) - sinf(att.phi)*sinf(att.theta)*sinf(att.psi))*Thrust_0; //(-9.81f);
		g22 = (cosf(att.phi)*cosf(att.theta)*sinf(att.psi))*Thrust_0; //(-9.81f);
		g23 = (-sinf(att.phi)*cosf(att.psi) + cosf(att.phi)*sinf(att.theta)*sinf(att.psi));
		g31 = (-sinf(att.phi)*cosf(att.theta))*Thrust_0; //(-9.81f);
		g32 = (-cosf(att.phi)*sinf(att.theta))*Thrust_0; //(-9.81f);
		g33 = (cosf(att.phi)*cosf(att.theta));
	}

	//Maybe add the found forumula for thrust here using the pwm and battery voltage data, such that it is not fixed around hover thrust

	// Next four blocks of the code are to compute the Moore-Penrose inverse of the G matrix
	// (G'*G)
	float a11 = g11*g11 + g21*g21 + g31*g31;
	float a12 = g11*g12 + g21*g22 + g31*g32;
	float a13 = g11*g13 + g21*g23 + g31*g33;
	float a21 = g12*g11 + g22*g21 + g32*g31;
	float a22 = g12*g12 + g22*g22 + g32*g32;
	float a23 = g12*g13 + g22*g23 + g32*g33;
	float a32 = g13*g12 + g23*g22 + g33*g32;
	float a33 = g13*g13 + g23*g23 + g33*g33;

	// Determinant of (G'*G)
	float detG = (a11*a22*a33 + a12*a23*a31 + a21*a32*a13) - (a13*a22*a31 + a11*a32*a23 + a12*a21*a33); 

	// Inverse of (G'*G)
	float a11_inv = (a22*a33 - a23*a32)/detG;
	float a12_inv = (a13*a32 - a12*a33)/detG;
	float a13_inv = (a12*a23 - a13*a22)/detG;
	float a21_inv = (a23*a31 - a21*a33)/detG;
	float a22_inv = (a11*a33 - a13*a31)/detG;
	float a23_inv = (a13*a21 - a11*a23)/detG;
	float a31_inv = (a21*a32 - a22*a31)/detG;
	float a32_inv = (a12*a31 - a11*a32)/detG;
	float a33_inv = (a11*a22 - a12*a21)/detG; 

	// G_inv = (G'*G)_inv*G'
	float g11_inv = a11_inv*g11 + a12_inv*g12 + a13_inv*g13;
	float g12_inv = a11_inv*g21 + a12_inv*g22 + a13_inv*g23;
	float g13_inv = a11_inv*g31 + a12_inv*g32 + a13_inv*g33;
	float g21_inv = a21_inv*g11 + a22_inv*g12 + a23_inv*g13;
	float g22_inv = a21_inv*g21 + a22_inv*g22 + a23_inv*g23;
	float g23_inv = a21_inv*g31 + a22_inv*g32 + a23_inv*g33;
	float g31_inv = a31_inv*g11 + a32_inv*g12 + a33_inv*g13;
	float g32_inv = a31_inv*g21 + a32_inv*g22 + a33_inv*g23;
	float g33_inv = a31_inv*g31 + a32_inv*g32 + a33_inv*g33;

	// Lin. accel. error multiplied  G^(-1) matrix (T_tilde negated because motor accepts only positiv commands, angles are in rad)
	if (ceiling_dist <= 0.1f){
		indiOuter.phi_tilde = (g11_inv*(indiOuter.linear_accel_err.x  - g14*K_ctr_eff*velS_z) + g12_inv*(indiOuter.linear_accel_err.y - g24*K_ctr_eff*velS_z) + g13_inv*(indiOuter.linear_accel_err.z - g34*K_ctr_eff*velS_z));
		indiOuter.theta_tilde = (g21_inv*(indiOuter.linear_accel_err.x  - g14*K_ctr_eff*velS_z) + g22_inv*(indiOuter.linear_accel_err.y - g24*K_ctr_eff*velS_z) + g23_inv*(indiOuter.linear_accel_err.z - g34*K_ctr_eff*velS_z));
		indiOuter.T_tilde = -(g31_inv*(indiOuter.linear_accel_err.x  - g14*K_ctr_eff*velS_z) + g32_inv*(indiOuter.linear_accel_err.y - g24*K_ctr_eff*velS_z) + g33_inv*(indiOuter.linear_accel_err.z - g34*K_ctr_eff*velS_z))/K_thr;
	}
	else{
		indiOuter.phi_tilde   = (g11_inv*indiOuter.linear_accel_err.x + g12_inv*indiOuter.linear_accel_err.y + g13_inv*indiOuter.linear_accel_err.z);
		indiOuter.theta_tilde = (g21_inv*indiOuter.linear_accel_err.x + g22_inv*indiOuter.linear_accel_err.y + g23_inv*indiOuter.linear_accel_err.z);
		indiOuter.T_tilde     = -(g31_inv*indiOuter.linear_accel_err.x + g32_inv*indiOuter.linear_accel_err.y + g33_inv*indiOuter.linear_accel_err.z)/K_thr; 	
	}

	// Filter thrust
	filter_thrust(indiOuter.thr, &indiOuter.T_incremented, &indiOuter.T_inner_f);

	// Pass thrust through the model of the actuator dynamics
	indiOuter.T_inner = indiOuter.T_inner + indiOuter.act_dyn_posINDI*(indiOuter.T_inner_f - indiOuter.T_inner); 

	// Compute trust that goes into the inner loop
	indiOuter.T_incremented = indiOuter.T_tilde + indiOuter.T_inner;

	// Compute commanded attitude to the inner INDI
	indiOuter.attitude_c.phi = indiOuter.attitude_f.phi + indiOuter.phi_tilde;
	indiOuter.attitude_c.theta = indiOuter.attitude_f.theta + indiOuter.theta_tilde;

	// Clamp commands
	indiOuter.T_incremented = clamp(indiOuter.T_incremented, MIN_THRUST, MAX_THRUST);
	indiOuter.attitude_c.phi = clamp(indiOuter.attitude_c.phi, -radians(pq_clamp), radians(pq_clamp));
	indiOuter.attitude_c.theta = clamp(indiOuter.attitude_c.theta, -radians(pq_clamp), radians(pq_clamp));

	// Reference values, which are passed to the inner loop INDI (attitude controller)
	refOuterINDI->x = radians(indiOuter.attitude_c.phi);
	refOuterINDI->y = radians(indiOuter.attitude_c.theta);
	refOuterINDI->z = indiOuter.T_incremented;
}

/**
 * Tuning settings for the gains of the INDI
 * controller for the position and velocity
 * of the Crazyflie in the X, Y and Z direction in the global
 * coordinate system.
 */
PARAM_GROUP_START(posCtrlIndi)

/**
 * @brief INDI position controller X proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_xi_x, &K_xi_x)
/**
 * @brief INDI position controller Y proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_xi_y, &K_xi_y)
/**
 * @brief INDI position controller Z proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_xi_z, &K_xi_z)

/**
 * @brief INDI velocity controller X proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_dxi_x, &K_dxi_x)
/**
 * @brief INDI velocity controller Y proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_dxi_y, &K_dxi_y)
/**
 * @brief INDI velocity controller Z proportional gain
 */
PARAM_ADD(PARAM_FLOAT, K_dxi_z, &K_dxi_z)

/**
 * @brief INDI Clamping value for the INDI roll and pitch command to inner loop [degrees]
 */
PARAM_ADD(PARAM_FLOAT, pq_clamping, &pq_clamp)

PARAM_GROUP_STOP(posCtrlIndi)

LOG_GROUP_START(posCtrlIndi)

/**
 * @brief INDI position reference input x [m]
 */
LOG_ADD(LOG_FLOAT, posRef_x, &positionRef.x)
/**
 * @brief INDI position reference input y [m]
 */
LOG_ADD(LOG_FLOAT, posRef_y, &positionRef.y)
/**
 * @brief INDI position reference input z [m]
 */
LOG_ADD(LOG_FLOAT, posRef_z, &positionRef.z)

/**
 * @brief INDI current velocity x [m/s]
 */
LOG_ADD(LOG_FLOAT, velS_x, &velS_x)
/**
 * @brief INDI current velocity y [m/s]
 */
LOG_ADD(LOG_FLOAT, velS_y, &velS_y)
/**
 * @brief INDI current velocity z [m/s]
 */
LOG_ADD(LOG_FLOAT, velS_z, &velS_z)

/**
 * @brief INDI velocity reference input x [m/s]
 */
LOG_ADD(LOG_FLOAT, velRef_x, &velocityRef.x)
/**
 * @brief INDI velocity reference input y [m/s]
 */
LOG_ADD(LOG_FLOAT, velRef_y, &velocityRef.y)
/**
 * @brief INDI velocity reference input z [m/s]
 */
LOG_ADD(LOG_FLOAT, velRef_z, &velocityRef.z)

/**
 * @brief INDI current attitude roll angle [rad]
 */
LOG_ADD(LOG_FLOAT, angS_roll, &indiOuter.attitude_s.phi)
/**
 * @brief INDI current attitude pitch angle [rad]
 */
LOG_ADD(LOG_FLOAT, angS_pitch, &indiOuter.attitude_s.theta)
/**
 * @brief INDI current attitude yaw angle [rad]
 */
LOG_ADD(LOG_FLOAT, angS_yaw, &indiOuter.attitude_s.psi)

/**
 * @brief INDI current attitude roll angle filtered (8 Hz low-pass) [rad]
 */
LOG_ADD(LOG_FLOAT, angF_roll, &indiOuter.attitude_f.phi)
/**
 * @brief INDI current attitude pitch angle filtered (8 Hz low-pass) [rad]
 */
LOG_ADD(LOG_FLOAT, angF_pitch, &indiOuter.attitude_f.theta)
/**
 * @brief INDI current attitude yaw angle filtered (8 Hz low-pass) [rad]
 */
LOG_ADD(LOG_FLOAT, angF_yaw, &indiOuter.attitude_f.psi)

/**
 * @brief INDI linear acceleration reference input x [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accRef_x, &indiOuter.linear_accel_ref.x)
/**
 * @brief INDI linear acceleration reference input y [m/s^2], NED frame 
 */
LOG_ADD(LOG_FLOAT, accRef_y, &indiOuter.linear_accel_ref.y)
/**
 * @brief INDI linear acceleration reference input z [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accRef_z, &indiOuter.linear_accel_ref.z)

/**
 * @brief INDI current linear acceleration measurement x [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accS_x, &indiOuter.linear_accel_s.x)
/**
 * @brief INDI current linear acceleration measurement z [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accS_y, &indiOuter.linear_accel_s.y)
/**
 * @brief INDI current linear acceleration measurement z [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accS_z, &indiOuter.linear_accel_s.z)

/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) x [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accF_x, &indiOuter.linear_accel_f.x)
/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) y [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accF_y, &indiOuter.linear_accel_f.y)
/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) z [m/s^2], body frame
 */
LOG_ADD(LOG_FLOAT, accF_z, &indiOuter.linear_accel_f.z)

/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) and rotated x [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accFT_x, &indiOuter.linear_accel_ft.x)
/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) and rotated y [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accFT_y, &indiOuter.linear_accel_ft.y)
/**
 * @brief INDI current linear acceleration measurement filtered (8 Hz low-pass) and rotated z [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accFT_z, &indiOuter.linear_accel_ft.z)

/**
 * @brief INDI linear acceleration error x [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accErr_x, &indiOuter.linear_accel_err.x)
/**
 * @brief INDI linear acceleration error y [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accErr_y, &indiOuter.linear_accel_err.y)
/**
 * @brief INDI linear acceleration error z [m/s^2], NED frame
 */
LOG_ADD(LOG_FLOAT, accErr_z, &indiOuter.linear_accel_err.z)

/**
 * @brief INDI roll angle command increment [rad]
 */
LOG_ADD(LOG_FLOAT, phi_tilde, &indiOuter.phi_tilde)
/**
 * @brief INDI pitch angle command increment [rad]
 */
LOG_ADD(LOG_FLOAT, theta_tilde, &indiOuter.theta_tilde)
/**
 * @brief INDI thrust command increment [motor units]
 */
LOG_ADD(LOG_FLOAT, T_tilde, &indiOuter.T_tilde)

/**
 * @brief INDI final previous thrust command, filtered with low-pass and passed through actuator dynamics
 */
LOG_ADD(LOG_FLOAT, T_inner, &indiOuter.T_inner)
/**
 * @brief INDI previous thrust command filtered (8 Hz low-pass) [motor units]
 */

LOG_ADD(LOG_FLOAT, T_inner_f, &indiOuter.T_inner_f)
/**
 * @brief INDI motor thrust command provided to inner loop [motor units]
 */
LOG_ADD(LOG_FLOAT, T_incremented, &indiOuter.T_incremented)
/**
 * @brief INDI roll angle command to inner loop [rad]
 */
LOG_ADD(LOG_FLOAT, cmd_phi, &indiOuter.attitude_c.phi)
/**
 * @brief INDI pitch angle command to inner loop [rad]
 */
LOG_ADD(LOG_FLOAT, cmd_theta, &indiOuter.attitude_c.theta)

LOG_ADD(LOG_FLOAT, T_0,  &Thrust_0)
LOG_ADD(LOG_FLOAT, T_0_N,  &Thrust_0_N)

LOG_ADD(LOG_FLOAT, motor_1_norm, &motor_norm_test[0])
LOG_ADD(LOG_FLOAT, motor_2_norm, &motor_norm_test[1])
LOG_ADD(LOG_FLOAT, motor_3_norm, &motor_norm_test[2])
LOG_ADD(LOG_FLOAT, motor_4_norm, &motor_norm_test[3])

LOG_ADD(LOG_UINT32, motor_1, &motorPower.m1)
LOG_ADD(LOG_UINT32, motor_2, &motorPower.m2)
LOG_ADD(LOG_UINT32, motor_3, &motorPower.m3)
LOG_ADD(LOG_UINT32, motor_4, &motorPower.m4)

LOG_ADD(LOG_FLOAT, ceil_dist, &ceiling_dist)
LOG_ADD(LOG_FLOAT, T_RBF_Ratio, &T_Ratio_RBF)
LOG_ADD(LOG_FLOAT, T_MEAS_Ratio, &T_RatioThrust_0_meas)

LOG_ADD(LOG_FLOAT, Weight_0, &RBF_t.weights_s[0])
LOG_ADD(LOG_FLOAT, Weight_1, &RBF_t.weights_s[1])
LOG_ADD(LOG_FLOAT, Weight_2, &RBF_t.weights_s[2])
LOG_ADD(LOG_FLOAT, Weight_3, &RBF_t.weights_s[3])
LOG_ADD(LOG_FLOAT, Weight_4, &RBF_t.weights_s[4])
LOG_ADD(LOG_FLOAT, Weight_5, &RBF_t.weights_s[5])
LOG_ADD(LOG_FLOAT, Weight_6, &RBF_t.weights_s[6])

LOG_ADD(LOG_FLOAT, ctr_eff_g14, &g14)
LOG_ADD(LOG_FLOAT, ctr_eff_g24, &g24)
LOG_ADD(LOG_FLOAT, ctr_eff_g34, &g34)



LOG_GROUP_STOP(posCtrlIndi)
