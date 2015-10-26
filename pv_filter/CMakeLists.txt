
SET(PLUGIN_NAME TimeAnnotate)

ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 
  SERVER_MANAGER_XML
    vtkDateTimeAnnotation.xml
  SERVER_MANAGER_SOURCES
    vtkDateTimeAnnotation.cxx
)
