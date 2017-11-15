/**
 * @file		collsisions.h
 * @brief		MCC Colisional module.
 * @author		Steffen Mattias Brask <steffemb@fys.uio.no>
 *
 * Functions to apply Monte Carlo Collisions between particles
 */

#ifndef COLLISIONS_H
#define COLLISIONS_H


/**
 * @brief Run designed for testing MCC
 * @param 	ini
 *
 *  Run designed for testing MCC
 *
 */
void mccTestMode(dictionary *ini);
funPtr mccTestMode_set(dictionary *ini);


/**
 * @brief brief description
 * @param one 		one description
 * @param two		two description
 *
 * Looong description of mccTest
 */


 //delete BorisTest in final Version
 void BorisTestMode(dictionary *ini);
 funPtr BorisTestMode_set(dictionary *ini);



int mccTest(int one, int two);

#endif // COLLISIONS_H
