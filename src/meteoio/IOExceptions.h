/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __IOEXCEPTIONS_H__
#define __IOEXCEPTIONS_H__

#include <exception>
#include <string>
#include <stdlib.h>

#include <meteoio/exports.h>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

namespace mio {

/**
 * @class IOException
 * @brief The basic exception class adjusted for the needs of SLF software
 *
 * @author Thomas Egger
 */


class IOException : public std::exception {
	public:
		IOException(const std::string& message="IOException occured", const std::string& position="");
		~IOException() throw() {}
		virtual const char* what() const throw();

	protected:
	#if defined(__linux) && !defined(ANDROID) && !defined(__CYGWIN__)
		std::string resolveSymbols(char *symbols, const unsigned int& ii, bool& found_main) const;
	#endif
		std::string msg, full_output;
};

/**
 * @class FileNotFoundException
 * @brief thrown when a there is an unsuccessful attempt to locate a file
 *
 * @author Thomas Egger
 */
class MIO_API FileNotFoundException : public IOException {
	public:
		FileNotFoundException(const std::string& filename="",
		                      const std::string& position="") : IOException("FileNotFoundException: " + filename, position){}
};

/**
 * @class FileAccessException
 * @brief thrown when a there are insufficient rights to access the file in a certain way (e.g. read, write)
 *
 * @author Thomas Egger
 */
class MIO_API FileAccessException : public IOException {
	public:
		FileAccessException(const std::string& filename="",
		                    const std::string& position="") : IOException("FileAccessException: " + filename, position){}
};

/**
 * @class InvalidFileNameException
 * @brief thrown when a filename given is not valid (e.g. "..", "." or empty)
 *
 * @author Thomas Egger
 */
class MIO_API InvalidFileNameException : public IOException {
	public:
		InvalidFileNameException(const std::string& filename="",
		                         const std::string& position="") : IOException("InvalidFileNameException: " + filename, position){}
};

/**
 * @class InvalidFormatException
 * @brief thrown when parsed data does not reflect an expected format (e.g. premature end of a line, file)
 *
 * @author Thomas Egger
 */
class MIO_API InvalidFormatException : public IOException {
	public:
		InvalidFormatException(const std::string& message="",
		                       const std::string& position="") : IOException("InvalidFormatException: " + message, position){}
};

/**
 * @class IndexOutOfBoundsException
 * @brief thrown when an index is out of bounds
 *
 * @author Thomas Egger
 */
class MIO_API IndexOutOfBoundsException : public IOException {
	public:
		IndexOutOfBoundsException(const std::string& message="",
		                          const std::string& position="") : IOException("IndexOutOfBoundsException: " + message, position){}
};

/**
 * @class ConversionFailedException
 * @brief thrown when an unsuccessful to convert data types/classes is made (e.g. attempt to convert a literal into a number)
 *
 * @author Thomas Egger
 */
class MIO_API ConversionFailedException : public IOException {
	public:
		ConversionFailedException(const std::string& message="",
		                          const std::string& position="") : IOException("ConversionFailedException: " + message, position){}
};

/**
 * @class InvalidArgumentException
 * @brief thrown when encountered an unexpected function's argument (e.g. bad index, bad or missing parameter name, etc.)
 *
 * @author Florian Hof
 */
class MIO_API InvalidArgumentException : public IOException {
	public:
		InvalidArgumentException(const std::string& message="",
		                         const std::string& position="") : IOException("InvalidArgumentException: " + message, position){}
};

/**
 * @class UnknownValueException
 * @brief thrown when encountered an unexpected value (e.g. unknown name or key)
 *
 * @author Florian Hof
 */
class MIO_API UnknownValueException : public IOException {
	public:
		UnknownValueException(const std::string& message="",
		                      const std::string& position="") : IOException("UnknownValueException: " + message, position){}
};

/**
 * @class NoAvailableDataException
 * @brief thrown when no data is available
 *
 * @author Florian Hof
 */
class MIO_API NoAvailableDataException : public IOException
{
	public:
		NoAvailableDataException(const std::string& message="",
		                         const std::string& position="") : IOException("NoAvailableDataException: " + message, position){}
};
} //end namespace

// Define DEBUG an empty function for seq compilation
#ifndef DEBUG
#define DEBUG printdebug
inline void printdebug(...) {}
#endif

#endif /*__IOException_H__*/


