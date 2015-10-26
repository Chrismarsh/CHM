/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile: vtkJulianDayToTextConvertor.h,v $

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDateTimeAnnotation
// .SECTION Description
// This filter can be attached to any filter/source/reader that supports time.
// vtkDateTimeAnnotation will generate a 1x1 vtkTable with the string
// for the data time using the format specified.
// The input to this filter is optional. If no input is specified, it will show
// produce request time in the output.

#pragma once

#include "vtkPVVTKExtensionsDefaultModule.h" //needed for exports
#include "vtkTableAlgorithm.h"

#include <boost/date_time/posix_time/posix_time.hpp>
using namespace boost::posix_time;

class VTKPVVTKEXTENSIONSDEFAULT_EXPORT vtkDateTimeAnnotation : public vtkTableAlgorithm
{
public:
  static vtkDateTimeAnnotation * New();
    vtkTypeMacro(vtkDateTimeAnnotation, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get/Set the format in which the to display the
  // input update time. Use printf formatting.
  // Default is "Time: %f".
  vtkSetStringMacro(Format);
  vtkGetStringMacro(Format);

  // Description:
  // Apply a translation to the time
  vtkSetMacro(Shift, double);
  vtkGetMacro(Shift, double);

  // Description:
  // Apply a scale to the time.
  vtkSetMacro(Scale, double);
  vtkGetMacro(Scale, double);

  // Description:
  // Convert the date/time from Gregorian numbers to Julian
  vtkSetMacro(ConvertFromGregorian, int);
  vtkGetMacro(ConvertFromGregorian, int);
  vtkBooleanMacro(ConvertFromGregorian, int);

// BTX
protected:
   vtkDateTimeAnnotation();
  ~vtkDateTimeAnnotation();

  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector);

  char  *Format;
  double Shift;
  double Scale;
  int    ConvertFromGregorian;

private:
  vtkDateTimeAnnotation(const vtkDateTimeAnnotation &); // Not implemented
  void operator=(const vtkDateTimeAnnotation &); // Not implemented
//ETX
};

