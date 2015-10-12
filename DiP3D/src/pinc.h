/**
 * @file		pinc.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC main routine.
 * @date		11.10.15
 *
 * PINC main header file. 
 * Holding the central function declarations comprising the PINC program.
 */

#ifndef PINC_H
#define PINC_H

/******************************************************************************
 * DEFINED IN INPUT.C
 *****************************************************************************/

/**
 * @brief	Parse PINC's input argument and input file
 * @param	argc	Argument count
 * @param	argv	Argument vector
 * @return	void 
 *
 * This function performs a sanity check of argc and argv and reads the
 * specified input file. It performs the necessary sanity checks of its
 * data and stores the values. It also computes derived values.
 * In case of failure it prints an error to stderr and terminates PINC.
 */
void parse_input(int argc, char *argv[]);

/******************************************************************************
 * DEFINED IN AUX.C
 *****************************************************************************/

/**
 * @brief	Terminates PINC due to error.
 * @param	format	Error message with specification of how to interpret data.
 * @param	...		Data to be interpreted in message.
 * @return	void 
 *
 * Similar syntax to printf(). Appends \n automatically at end of line.
 */
void pincerror(const char* restrict format,...);

#endif // PINC_H
