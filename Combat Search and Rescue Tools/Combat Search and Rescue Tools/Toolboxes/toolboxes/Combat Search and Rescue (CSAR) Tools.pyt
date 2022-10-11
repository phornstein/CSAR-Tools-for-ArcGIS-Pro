"""
TITLE: Combat Search and Rescue (CSAR) Tools
DATE: September 2022
VERSION: 1.0
DEVELOPER: Chris Lee (clee@esri.com), Parker Hornstein (phornstein@esri.com)
REQUIREMENTS: ArcGIS Pro 2.9.x, Spatial Analyst License
LICENSE:

Copyright © 2022 Esri

All rights reserved under the copyright laws of the United States and applicable international laws, treaties, and conventions.
You may freely redistribute and use this sample code, with or without modification, provided you include the original copyright notice and use restrictions.

Disclaimer: THE SAMPLE CODE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ESRI OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
SUSTAINED BY YOU OR A THIRD PARTY, HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT ARISING IN ANY WAY OUT OF
THE USE OF THIS SAMPLE CODE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For additional information, contact:

Esri
Attn: Contracts and Legal Services Department
380 New York Street
Redlands, California, 92373-8100
USA
email: contracts@esri.com
"""


import os
import arcpy
from arcpy import env
from math import radians,sin,cos

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Combat Search and Rescue (CSAR) Tools"
        self.alias = "Combat Search and Rescue (CSAR) Tools"

        # List of tool classes associated with this toolbox
        self.tools = [DefineAOI, HastyHLZSuitability, DeliberateHLZSuitability, CreateOperationsGraphic, FishnetHLZAnalysis]


#---------- Define AOI Tool ----------#

class DefineAOI(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Define AOI"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        poi = arcpy.Parameter(
            displayName="9-Line Location:",
            name="poi",
            datatype="GPFeatureRecordSetLayer",
            parameterType="Optional",
            direction="Input"
        )

        aoi = arcpy.Parameter(
            displayName="Area of Operations:",
            name="aoi",
            datatype="GPFeatureRecordSetLayer",
            parameterType="Optional",
            direction="Input"
        )

        mgrs = arcpy.Parameter(
            displayName="9-Line Location MGRS:",
            name="mgrs",
            datatype="GPString",
            parameterType="Optional",
            direction="Input"
        )
        
        mgrs_buffer = arcpy.Parameter(
            displayName="Buffer Distance for MGRS in Meters:",
            name="mgrs_buffer",
            datatype="GPLinearUnit",
            parameterType="Optional",
            direction="Input"
        )

        op_name = arcpy.Parameter(
            displayName="Operation Name:",
            name="op_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )

        out_point = arcpy.Parameter(
            displayName="Point Featureclass:",
            name="out_point",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output"
        )

        out_aoi = arcpy.Parameter(
            displayName="AOI Featureclass:",
            name="out_aoi",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output"
        )

        poi.filter.list = ['Point']
        aoi.filter.list = ['Polygon']

        out_point.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\aoi_point.lyrx')
        out_aoi.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\aoi_polygon.lyrx')

        params = [poi,aoi,mgrs,mgrs_buffer,op_name,out_point,out_aoi]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True
    
    def _disableParam(self,parameter):
        parameter.value = None
        parameter.enabled = False

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if not parameters[0].value is None: #if point is drawn disable MGRS
            self._disableParam(parameters[2])
        else:
            parameters[2].enabled = True

        if not parameters[2].value is None: #if MGRS is entered disabled point drawing
            self._disableParam(parameters[0])
        else:
            parameters[0].enabled = True

        if not parameters[1].value is None: #if an AOI is drawn disable buffer
            self._disableParam(parameters[3])
        else:
            parameters[3].enabled = True

        if not parameters[3].value is None: #if buffer distance is enabled disable AOI draw
            self._disableParam(parameters[1])
        else:
            parameters[3].value = None
            parameters[3].enabled = True

        if not parameters[4].value is None and bool(parameters[5].value) and bool(parameters[6].value):
            op_name = parameters[4].valueAsText
            pnt_name = os.path.join(arcpy.env.workspace,arcpy.ValidateTableName("{}_point".format(op_name)))
            parameters[5].value = pnt_name

            aoi_name = os.path.join(arcpy.env.workspace,arcpy.ValidateTableName("{}_aoi".format(op_name)))
            parameters[6].value = aoi_name
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        poi = parameters[0].valueAsText
        aoi = parameters[1].valueAsText
        mgrs = parameters[2].valueAsText
        mgrs_buffer = parameters[3].value
        op_name = parameters[4].valueAsText
        out_point = parameters[5].valueAsText
        out_aoi = parameters[6].valueAsText
        sr = arcpy.SpatialReference(4326)

        #create point
        arcpy.SetProgressorLabel("Creating Point Featureclass...")
        if not poi is None:
            arcpy.CopyFeatures_management(poi,out_point)
            arcpy.Delete_management(poi) #delete the temp layer mostly cause i was annoyed with it still being there
        else:
            geom = arcpy.FromCoordString(mgrs,"MGRS")
            geom = geom.firstPoint

            fc = arcpy.CreateFeatureclass_management(os.path.dirname(out_point),os.path.basename(out_point),"POINT",spatial_reference=sr)
            with arcpy.da.InsertCursor(fc,['SHAPE@']) as cursor:
                cursor.insertRow([geom])

        #create AOI
        arcpy.SetProgressorLabel("Creating Polygon Featureclass...")
        if not aoi is None:
            arcpy.CopyFeatures_management(aoi,out_aoi)
            arcpy.Delete_management(aoi) #delete the temp layer mostly cause i was annoyed with it still being there
        else:
            res = arcpy.Buffer_analysis(out_point,arcpy.Geometry(),mgrs_buffer) #write the buffer directly to a polygon geometry
            extent_poly = res[0].extent.polygon #get the polygon extent of the buffer
            fc = arcpy.CreateFeatureclass_management(os.path.dirname(out_aoi),os.path.basename(out_aoi),"POLYGON",spatial_reference=sr)
            with arcpy.da.InsertCursor(fc,['SHAPE@']) as cursor:
                cursor.insertRow([extent_poly])

        return



#---------- Hasty HLZ Suitability Tool ----------#

class HastyHLZSuitability(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Hasty HLZ Suitability"
        self.description = "This tool will generate an HLZ Suitability layer for a given MEDEVAC location based on elevation, landcover, and vertical obstruction (DVOF) data."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        # AOI
        aoi = arcpy.Parameter(
            displayName = "AOI Polygon",
            name = "aoi",
            datatype = "GPFeatureRecordSetLayer",
            parameterType = "Required",
            direction = "Input")
        aoi.filter.list = ["Polygon"]
        
        # ELEVATION LAYER
        elevationParam = arcpy.Parameter(
            displayName = "Elevation Surface",
            name = "elevation",
            datatype = "GPRasterLayer",
            parameterType = "Optional",
            direction = "Input")

        # SLOPE LAYER
        slopeParam = arcpy.Parameter(
            displayName = "Slope Surface",
            name = "slopeParam",
            datatype = "GPRasterLayer",
            parameterType = "Optional",
            direction = "Input")

        # MAX SLOPE
        maxSlope = arcpy.Parameter(
            displayName="Maximum Allowable Slope:",
            name="maxSlope",
            datatype="GPLong",
            parameterType="Optional",
            direction="Input")
        maxSlope.value = 15

        # LANDCOVER LAYER
        landcoverParam = arcpy.Parameter(
            displayName = "Landcover Surface",
            name = "landcover",
            datatype = "GPRasterLayer",
            parameterType = "Required",
            direction = "Input")

        # LANDCOVER TYPE
        landcoverTypeParam = arcpy.Parameter(
            displayName = "Landcover Type",
            name = "landcoverType",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        # set up parameter's value list for different landcover types
        landcoverTypeParam.filter.type = "ValueList"
        landcoverTypeParam.filter.list = ["NLCD", "GeoCover", "VISNAV"]

        # BUILDING OBSTRUCTIONS LAYER
        buildingObstructionsParam = arcpy.Parameter(
            displayName = "Horizontal Obstruction Layers",
            name = "buildingObstructionsParam",
            datatype = "GPFeatureLayer",
            parameterType = "Optional",
            direction = "Input",
            multiValue=True)
        buildingObstructionsParam.filter.list = ['Polygon']    

        # VERTICAL OBSTRUCTIONS LAYER
        verticalObstructionsParam = arcpy.Parameter(
            displayName = "Vertical Obstructions Layer",
            name = "verticalObstructionsParam",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        verticalObstructionsParam.filter.list = ['Point']

        # VERTICAL OBSTRUCTIONS SOURCE
        verticalObstructionBufferField = arcpy.Parameter(
            displayName = "Vertical Obstructions Buffer Field",
            name = "verticalObstructionBufferField",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

         # HLZ SUITABILITY OUTPUT
        hlzSuitabilityOutputParam = arcpy.Parameter(
            displayName = "HLZ Polygon Candidates",
            name = "hlzSuitabilityOutput",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Output")

         # HLZ SUITABILITY RASTER OUTPUT
        hlzSuitabilityRasterOutputParam = arcpy.Parameter(
            displayName = "HLZ Suitability Raster Output",
            name = "hlzSuitabilityRasterOutputParam",
            datatype = "DERasterDataset",
            parameterType = "Optional",
            direction = "Output")

        params = [aoi, 
                elevationParam,
                slopeParam,
                maxSlope,
                landcoverParam, 
                landcoverTypeParam, 
                buildingObstructionsParam,
                verticalObstructionsParam, 
                verticalObstructionBufferField, 
                hlzSuitabilityOutputParam,
                hlzSuitabilityRasterOutputParam]

        # set symbology for HLZ Suitability
        hlzSuitabilityOutputParam.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\hlz_polys.lyrx')
        hlzSuitabilityRasterOutputParam.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\hlz_raster.lyrx')

        return params

    def isLicensed(self):
        """Allow the tool to execute, only if the ArcGIS Spatial Analyst extension is available."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed

        return True  # tool can be executed

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if not parameters[1].value is None:
            parameters[2].parameterType = 'Optional'
            parameters[2].enabled = False
        elif not parameters[2].value is None:
            parameters[1].parameterType = 'Optional'
            parameters[1].enabled = False
        else:
            parameters[1].enabled = True
            parameters[2].enabled = True

        if not parameters[7].value is None:
            flds = [f.name for f in arcpy.ListFields(parameters[7].valueAsText) if f.type in ['String','Integer','Double']]
            #update field input
            parameters[8].filter.list = flds

        if not parameters[0].value is None and (parameters[9].value is None and parameters[9].value is None):
            parameters[9].value = parameters[0].valueAsText + '_HLZ_Candidates'
            parameters[10].value = parameters[0].valueAsText + '_HLZ_Raster'
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        #check for elevation data
        if not bool(parameters[1].value) and not bool(parameters[2].value):
            arcpy.AddError("Either elevation or slope data is required....")
            exit(0)

        aoi = parameters[0].valueAsText
        elevationRaster = parameters[1].valueAsText if not parameters[1].value is None else None
        slopeRaster = parameters[2].value if not parameters[2].value is None else None
        maxSlope = parameters[3].valueAsText
        landcoverRaster = parameters[4].valueAsText
        landcoverType = parameters[5].valueAsText
        buildingLayer = parameters[6].values
        if buildingLayer is not None:
            buildingLayer = ';'.join([b.dataSource for b in buildingLayer])
        vofLayer = parameters[7].valueAsText
        vofBufferField = parameters[8].valueAsText
        output = parameters[9].valueAsText
        rasterOutput = parameters[10].valueAsText
      
        with arcpy.EnvManager(extent=arcpy.Describe(aoi).extent):
            arcpy.SetProgressorLabel("Buffering vertical obstructions...")
            vofBuffer = r'in_memory\vofBuffer'
            arcpy.analysis.Buffer(vofLayer,vofBuffer,vofBufferField,"FULL","ROUND","ALL",None,"PLANAR")

            if not buildingLayer is None:
                arcpy.SetProgressorLabel("Merging building obstructions...")
                vof = r'in_memory\vof'
                arcpy.Merge_management('{};{}'.format(buildingLayer,vofBuffer),vof)
            else:
                vof = vofBuffer

            vofFinalBuffer = r'in_memory\finalVofBuffer'
            arcpy.Buffer_analysis(vof,vofFinalBuffer,"12.5 Meters") #buffer all obstructions by 12.5m... one half the minimum distance for a helicopter land

            arcpy.management.CalculateField(vofFinalBuffer,"value","1","PYTHON3",'',"SHORT","NO_ENFORCE_DOMAINS")

            arcpy.SetProgressorLabel("Converting vertical obstructions to raster...")
            vofRaster = r'in_memory\vofRaster'
            arcpy.conversion.PolygonToRaster(vofFinalBuffer, "value", vofRaster, "CELL_CENTER", "NONE")

            arcpy.SetProgressorLabel("Reclassifying vertical obstructions")
            vofReclass = r'in_memory\vofReclass'
            result = arcpy.sa.Reclassify(vofRaster, "Value", "0 1;0 1 0;NODATA 1", "NODATA")
            result.save(vofReclass)

            if elevationRaster is not None:
                arcpy.SetProgressorLabel("Creating slope from elevation...")
                slopeRaster = r'in_memory\slopeRaster'
                result = arcpy.sa.Slope(elevationRaster, "DEGREE")
                result.save(slopeRaster)

            arcpy.SetProgressorLabel("Reclassifying slope...")
            arcpy.AddMessage("Slope of 0 -> {} degrees is acceptable".format(maxSlope))
            slopeReclass = r'in_memory\slopeRaster'
            remap = "0 {0} 1;{0} 90 0;NODATA 0".format(maxSlope)
            result = arcpy.sa.Reclassify(slopeRaster, "Value", remap, "NODATA")
            result.save(slopeReclass)

            arcpy.SetProgressorLabel("Reclassifying landcover...")
            landcoverReclass = r'in_memory\landcoverReclass'
            if landcoverType == "NLCD":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "11 0;21 0;23 0;24 0;31 1;41 0;42 0;43 0;51 1;52 1;71 1;72 1;73 1;74 1;81 1;82 1;90 0;95 0", "NODATA")
            elif landcoverType == "GeoCover":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "1 0;2 0;3 1;4 1;5 1;6 0;7 1;8 1;9 0;10 0;12 1;14 0;15 0", "NODATA")
            elif landcoverType == "VISNAV":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "1 0;2 0;3 1;4 1;5 1;6 0;7 1;8 1;9 0;10 0;11 0;12 1;14 1;15 1;16 1;17 1;18 1;19 1;20 0;21 0;27 1;28 0; 29 1", "NODATA")
            else:
                arcpy.AddError("Landcover reclass failed")
                exit(0)
            result.save(landcoverReclass)

            arcpy.SetProgressorLabel("Calculating HLZ...")
            vofXSlope = r'in_memory\vofXSlope'
            result = arcpy.sa.Times(vofReclass,slopeReclass)
            result.save(vofXSlope)

            combined = r'in_memory\combined'
            result = arcpy.sa.Times(vofXSlope,landcoverReclass)
            result.save(combined)

            finalRaster = r'in_memory\finalRaster'
            result = arcpy.sa.Reclassify(combined, "Value", "0 1;0 1 0;NODATA 1", "NODATA")
            result.save(finalRaster)
            if not rasterOutput is None:
                arcpy.SetProgressorLabel("Saving HLZ Raster...")
                result.save(rasterOutput)
            
            arcpy.SetProgressorLabel("Creating HLZ Candidate Polygons...")
            arcpy.AddMessage("HLZ candidate sizes are based on MCRP 2-10B.5")
            polys = r'in_memory\polys'
            arcpy.conversion.RasterToPolygon(finalRaster,polys, "NO_SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)

            selection = arcpy.management.SelectLayerByAttribute(polys, "NEW_SELECTION", "gridcode = 0", None)
            arcpy.CopyFeatures_management(selection,output)

            arcpy.management.CalculateGeometryAttributes(output, "area AREA_GEODESIC", '', "SQUARE_METERS", None, "SAME_AS_INPUT")           
            arcpy.management.CalculateField(output, "lz_size", "foo(!area!)", "PYTHON3", """def foo(area):
    if area < 625:
        return "Too Small"
    if 625 < area < 1225:
        return "Small LP"
    if 1225 < area < 2500:
        return "Medium LP"
    else:
        return "Large LP" """, "TEXT", "NO_ENFORCE_DOMAINS")

            arcpy.management.DeleteField(output, "gridcode;ORIG_FID;MBG_APodX1;MBG_APodY1;MBG_APodX2;MBG_APodY2", "DELETE_FIELDS")

        return



#---------- Deliberate HLZ Suitability Tool ----------#

class DeliberateHLZSuitability(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Deliberate HLZ Suitability"
        self.description = "This tool will identify three potential HLZs based off of the suitability, airframe landing zone requirements, and proximity to the MEDEVAC location."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # AOI
        aoi = arcpy.Parameter(
            displayName = "AOI Polygon",
            name = "aoi",
            datatype = "GPFeatureRecordSetLayer",
            parameterType = "Required",
            direction = "Input")
        aoi.filter.list = ["Polygon"]

        # AIRFRAME
        airframeParam = arcpy.Parameter(
            displayName = "Airframe",
            name = "airframe",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        airframe_table = os.path.join(os.path.dirname(__file__),r'HLZ_Toolbox_Development.gdb\Aircraft_Specifications') #change when this becomes a template
        field = 'Model'
        
        airframeParam.filter.type = "ValueList"
        airframeParam.filter.list = sorted([row[0] for row in arcpy.da.SearchCursor(airframe_table, field)])

        # ELEVATION
        elevationParam = arcpy.Parameter(
            displayName = "Elevation Surface",
            name = "elevation",
            datatype = "GPRasterLayer",
            parameterType = "Optional",
            direction = "Input")

        # SLOPE
        slopeParam = arcpy.Parameter(
            displayName = "Slope Surface",
            name = "slopeParam",
            datatype = "GPRasterLayer",
            parameterType = "Optional",
            direction = "Input")

        # LANDCOVER
        landcoverParam = arcpy.Parameter(
            displayName = "Landcover Surface",
            name = "landcover",
            datatype = "GPRasterLayer",
            parameterType = "Required",
            direction = "Input")

        # LANDCOVER TYPE
        landcoverTypeParam = arcpy.Parameter(
            displayName = "Landcover Type",
            name = "landcoverType",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        # set up parameter's value list for different landcover types
        landcoverTypeParam.filter.type = "ValueList"
        landcoverTypeParam.filter.list = ["NLCD", "GeoCover", "VISNAV"]

        # BUILDING/HORIZONTAL OBSTRUCTIONS
        buildingObstructionsParam = arcpy.Parameter(
            displayName = "Horizontal Obstruction Layers",
            name = "buildingObstructionsParam",
            datatype = "GPFeatureLayer",
            parameterType = "Optional",
            direction = "Input",
            multiValue=True)
        buildingObstructionsParam.filter.list = ['Polygon']    

        # VERTICAL OBSTRUCTIONS
        verticalObstructionsParam = arcpy.Parameter(
            displayName = "Vertical Obstructions Layer",
            name = "verticalObstructionsParam",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        verticalObstructionsParam.filter.list = ['Point']

        # VERTICAL OBSTRUCTIONS SOURCE
        verticalObstructionBufferField = arcpy.Parameter(
            displayName = "Vertical Obstructions Buffer Field",
            name = "verticalObstructionBufferField",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")

        # REPLACEMENT FOR HLZ SUITABILITY POLYGON OUTPUT
        hlzPolygonOutputParam = arcpy.Parameter(
            displayName = "HLZ Polygons",
            name = "hlzPolygonsOutput",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Output")

        # HLZ POINTS OUTPUT
        hlzPointsOutputParam = arcpy.Parameter(
            displayName = "HLZ Points",
            name = "hlzPointsOutput",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Output")

        # return list of parameters
        params = [aoi,
                airframeParam, 
                elevationParam,
                slopeParam,
                landcoverParam, 
                landcoverTypeParam, 
                buildingObstructionsParam,
                verticalObstructionsParam, 
                verticalObstructionBufferField, 
                hlzPolygonOutputParam,
                hlzPointsOutputParam]

        # set symbology for polygon and point outputs
        hlzPolygonOutputParam.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\deliberate_hlz_polygons.lyrx')
        hlzPointsOutputParam.symbology = os.path.join(os.path.dirname(__file__),r'layer_templates\deliberate_hlz_points.lyrx')

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # if elevation layer is provided, disable slope parameter; and vice versa
        if not parameters[2].value is None:
            parameters[3].parameterType = 'Optional'
            parameters[3].enabled = False
        elif not parameters[3].value is None:
            parameters[2].parameterType = 'Optional'
            parameters[2].enabled = False
        else:
            parameters[2].enabled = True
            parameters[3].enabled = True

        # when vertical obstruction layer provided, populate list of fields for that layer
        if not parameters[7].value is None:
            flds = [f.name for f in arcpy.ListFields(parameters[7].valueAsText) if f.type in ['String','Integer','Double']]
            #update field input
            parameters[8].filter.list = flds

        # adjust names of outputs
        if not parameters[0].value is None and (parameters[9].value is None and parameters[9].value is None):
            parameters[9].value = parameters[0].valueAsText + '_HLZ_Polygons'
            parameters[10].value = parameters[0].valueAsText + '_HLZ_Points'
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # check for elevation data, exit if not provided
        if not bool(parameters[1].value) and not bool(parameters[2].value):
            arcpy.AddError("Either elevation or slope data is required....")
            exit(0)
        
        # get parameters
        aoi = parameters[0].valueAsText
        airframe = parameters[1].valueAsText
        elevationRaster = parameters[2].valueAsText if not parameters[2].value is None else None
        slopeRaster = parameters[3].value if not parameters[3].value is None else None
        landcoverRaster = parameters[4].valueAsText
        landcoverType = parameters[5].valueAsText
        buildingLayer = parameters[6].values
        if buildingLayer is not None:
            buildingLayer = ';'.join([b.dataSource for b in buildingLayer])
        vofLayer = parameters[7].valueAsText
        vofBufferField = parameters[8].valueAsText
        hlzPolygonOutput = parameters[9].valueAsText
        hlzPointOutput = parameters[10].valueAsText

        # get airframe specifications
        airframe_table = os.path.join(os.path.dirname(__file__),r'HLZ_Toolbox_Development.gdb\Aircraft_Specifications')
        fields = ['Model', 'TDP_Diameter_Meters', 'Ground_Slope_Max_Percent']
        with arcpy.da.SearchCursor(airframe_table, fields , "Model = '{}'".format(airframe)) as cursor:
            for row in cursor:
                td_area = row[1]**2
                td_maxslope = row[2]
                break
            else:
                arcpy.AddError("Airframe specifications not found....")

        # begin processing
        arcpy.AddMessage("Starting 'Deliberate HLZ Suitability' for {}...".format(airframe))

        # set the processing extent to the extent of the 'aoi' layer
        with arcpy.EnvManager(extent=arcpy.Describe(aoi).extent):

            # buffer vertical obstructions
            arcpy.SetProgressorLabel("Buffering vertical obstructions...")
            vofBuffer = r'in_memory\vofBuffer'
            arcpy.analysis.Buffer(vofLayer,vofBuffer,vofBufferField,"FULL","ROUND","ALL",None,"PLANAR")

            # merge with horizontal obstructions (i.e. roads) if provided
            if not buildingLayer is None:
                arcpy.SetProgressorLabel("Merging building obstructions...")
                vof = r'in_memory\vof'
                arcpy.Merge_management('{};{}'.format(buildingLayer,vofBuffer),vof)
            else:
                vof = vofBuffer

            vofFinalBuffer = r'in_memory\finalVofBuffer'
            arcpy.Buffer_analysis(vof,vofFinalBuffer,"12.5 Meters") #buffer all obstructions by 12.5m... one half the minimum distance for a helicopter land

            # convert obstruction polygons to raster
            arcpy.management.CalculateField(vofFinalBuffer,"value","1","PYTHON3",'',"SHORT","NO_ENFORCE_DOMAINS")
            arcpy.SetProgressorLabel("Converting vertical obstructions to raster...")
            vofRaster = r'in_memory\vofRaster'
            arcpy.conversion.PolygonToRaster(vofFinalBuffer, "value", vofRaster, "CELL_CENTER", "NONE")

            # reclassify raster
            arcpy.SetProgressorLabel("Reclassifying vertical obstructions")
            vofReclass = r'in_memory\vofReclass'
            result = arcpy.sa.Reclassify(vofRaster, "Value", "0 1;0 1 0;NODATA 1", "NODATA")
            result.save(vofReclass)

            # if elavation is provided instead of slope, calculate slope
            if elevationRaster is not None:
                arcpy.SetProgressorLabel("Creating slope from elevation...")
                slopeRaster = r'in_memory\slopeRaster'
                result = arcpy.sa.Slope(elevationRaster, "DEGREE")
                result.save(slopeRaster)

            # reclassify slope
            arcpy.SetProgressorLabel("Reclassifying slope...")
            arcpy.AddMessage("Slope of 0 -> {} degrees is acceptable".format(td_maxslope))
            slopeReclass = r'in_memory\slopeRaster'
            remap = "0 {0} 1;{0} 90 0;NODATA 0".format(td_maxslope)
            result = arcpy.sa.Reclassify(slopeRaster, "Value", remap, "NODATA")
            result.save(slopeReclass)

            # reclassify landcover
            arcpy.SetProgressorLabel("Reclassifying landcover...")
            landcoverReclass = r'in_memory\landcoverReclass'
            if landcoverType == "NLCD":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "11 0;21 0;23 0;24 0;31 1;41 0;42 0;43 0;51 1;52 1;71 1;72 1;73 1;74 1;81 1;82 1;90 0;95 0", "NODATA")
            elif landcoverType == "GeoCover":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "1 0;2 0;3 1;4 1;5 1;6 0;7 1;8 1;9 0;10 0;12 1;14 0;15 0", "NODATA")
            elif landcoverType == "VISNAV":
                result = arcpy.sa.Reclassify(landcoverRaster, "VALUE", "1 0;2 0;3 1;4 1;5 1;6 0;7 1;8 1;9 0;10 0;11 0;12 1;14 1;15 1;16 1;17 1;18 1;19 1;20 0;21 0;27 1;28 0; 29 1", "NODATA")
            else:
                arcpy.AddError("Landcover reclass failed")
                exit(0)
            result.save(landcoverReclass)

            # calculate "obstructions" * "slope" raster
            arcpy.SetProgressorLabel("Calculating HLZ...")
            vofXSlope = r'in_memory\vofXSlope'
            result = arcpy.sa.Times(vofReclass,slopeReclass)
            result.save(vofXSlope)

            # calculate "combined" raster (obstructions * slope * landcover)
            combined = r'in_memory\combined'
            result = arcpy.sa.Times(vofXSlope,landcoverReclass)
            result.save(combined)

            # reclassify "combined" raster to get final raster output
            finalRaster = r'in_memory\finalRaster'
            result = arcpy.sa.Reclassify(combined, "Value", "0 1;0 1 0;NODATA 1", "NODATA")
            result.save(finalRaster)
            
            # create polygons
            arcpy.SetProgressorLabel("Creating HLZ Polygons...")
            arcpy.AddMessage("HLZ candidate sizes are based on MCRP 2-10B.5")
            polys = r'in_memory\polys'
            arcpy.conversion.RasterToPolygon(finalRaster,polys, "NO_SIMPLIFY", "Value", "SINGLE_OUTER_PART", None)

            hlz_polygons = r'in_memory\hlz_polygons'
            selection = arcpy.management.SelectLayerByAttribute(polys, "NEW_SELECTION", "gridcode = 0", None)
            arcpy.CopyFeatures_management(selection,hlz_polygons)

            # calculate geodesic area
            arcpy.management.CalculateGeometryAttributes(hlz_polygons, "area AREA_GEODESIC", '', "SQUARE_METERS", None, "SAME_AS_INPUT")           

            # delete unnecessary fields
            arcpy.management.DeleteField(hlz_polygons, "gridcode;ORIG_FID;MBG_APodX1;MBG_APodY1;MBG_APodX2;MBG_APodY2", "DELETE_FIELDS")

            # select polygons that meet area requirements for airframe
            arcpy.AddMessage("Searching for HLZs that are at least {} square meters.".format(td_area))
            area_selection = arcpy.management.SelectLayerByAttribute(hlz_polygons, "NEW_SELECTION", "area >= {}".format(td_area), None) #TODO area calculation to reflect the airframe's specs
            
            # generate Near Table to find three closest polygons to the MEDEVAC location (center of aoi) - use .5 distance of AOI polygon
            nears_point = arcpy.management.FeatureToPoint(aoi, r"in_memory\nears_point")
            nears_table = r"in_memory\nears_table"
            arcpy.analysis.GenerateNearTable(nears_point, area_selection, nears_table, "1500 Meters", "NO_LOCATION", "NO_ANGLE", "ALL", 3, "PLANAR") #TODO get polygon length or width for search distance
            
            # select polygons whose OBJECTIDs are in the list of values from the 'NEAR_FID' field in the Near Table
            featureids_in_nears_table = [row[0] for row in arcpy.da.SearchCursor(nears_table, "NEAR_FID")]
            nearest_selection = arcpy.management.SelectLayerByAttribute(area_selection, "NEW_SELECTION", "OBJECTID IN {}".format(tuple(featureids_in_nears_table)), None)
            
            # save polygons
            arcpy.CopyFeatures_management(nearest_selection, hlzPolygonOutput)

            # generate centroid points in those polygons
            arcpy.management.FeatureToPoint(hlzPolygonOutput, hlzPointOutput, "INSIDE")
            arcpy.management.CalculateGeometryAttributes(hlzPointOutput, [["mgrs", "POINT_COORD_NOTATION"]], "", "", "", "MGRS")


        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return



#---------- Create Operations Graphic (Layout) Tool ----------#

class CreateOperationsGraphic(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Create Operations Graphic (Layout)"
        self.description = "Create operational graphic layers and create a layout."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        # MGRS Coordinate
        layout_template = arcpy.Parameter(
            displayName = "Layout Template",
            name = "layout_template",
            datatype = "DEFile",
            parameterType = "Required",
            direction = "Input")

        layer_test_extent = arcpy.Parameter(
            displayName = "Layer",
            name = "layer",
            datatype = "GPLayer",
            parameterType = "Required",
            direction = "Input")

        graphic_name = arcpy.Parameter(
            displayName = "Graphic Name",
            name = "graphic_name",
            datatype = "String",
            parameterType = "Optional",
            direction = "Input")
        graphic_name.value = "Operations Graphic"

        params = [layout_template, layer_test_extent, graphic_name]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        #if not arcpy.Describe(elevation).spatialReference == arcpy.Describe(landcover).spatialReference:

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
         # Current Project and Map
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        active_map = aprx.activeMap
        
        # Reference layout template
        layout_template = parameters[0].valueAsText

        # Create layout object
        aprx.importDocument(layout_template, reuse_existing_maps=True)
        layout = aprx.listLayouts("Layout")[0]
        layout.name = parameters[2].valueAsText
        
        # Change the 'active' map in the layout's mapframe element
        mapframe = layout.listElements("MAPFRAME_ELEMENT", "Map Frame")[0] # return the 'MAPFRAME_ELEMENT' object titled "Map Frame" (see all the layout's elements in the Drawing Order on the Content pane when in Layout view)
        mapframe.map = active_map

        # Set the extent of the mapframe element
        lyr = parameters[1].valueAsText
        mapframe.camera.setExtent(arcpy.Describe(lyr).extent)

        # Set the title of the graphic
        title = layout.listElements("TEXT_ELEMENT", "Title")[0] # return the 'TEXT_ELEMENT' object titled "Title" (see all the layout's elements in the Drawing Order on the Content pane when in Layout view)
        title.text = parameters[2].valueAsText

        # Export layout as pdf
        products_dir = os.path.join(os.path.dirname(aprx.filePath),r'products') # directory to store exported pdfs
        if not os.path.exists(products_dir): # boolean to see if "products" folder exists - if not then create one
            os.makedirs(products_dir)
        layout.exportToPDF(os.path.join(products_dir,layout.name + '.pdf')) # export the layout
        
        return



#---------- Fishnet HLZ Analysis Testing ----------#

class FishnetHLZAnalysis(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Fishnet HLZ Analysis (still in development)"
        self.description = "Fill this out later.."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""

        input_features = arcpy.Parameter(
            displayName = "Input",
            name = "input_features",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")

        output_features = arcpy.Parameter(
            displayName = "Output",
            name = "output_features",
            datatype = "DEFeatureClass",
            parameterType = "Required",
            direction = "Output")

        params = [input_features, output_features]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        #if not arcpy.Describe(elevation).spatialReference == arcpy.Describe(landcover).spatialReference:

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        # Start tool here
        infc=parameters[0].valueAsText
        outFC=parameters[1].valueAsText

        d=arcpy.Describe(infc)
        SR=d.spatialReference

        # get size of polygons based on airframe specifications
        W=80;L=80;A=0.99*W*L
        fnet="in_memory/fnet"
        erased="in_memory/fnet"
        
        # function for rotating the polygon
        def ShapeMake(pGon,angle):
            a=radians(angle)
            ARR=arcpy.Array()
            cX=cPoint.X;cY=cPoint.Y
            for part in pGon.boundary():
                ar=arcpy.Array()
                for p in part:
                    x,y=p.X-cX,p.Y-cY
                    xN=cos(a)*x+sin(a)*y
                    yN=-sin(a)*x+cos(a)*y
                    pN=arcpy.Point(xN+cX,yN+cY)
                    ar.add(pN)
                ARR.add(ar)
            pgonRotated=arcpy.Polygon(ARR,SR)
            return pgonRotated

        # function for creating fishnet and counting complete polygons
        def fnetMake():
            FNET=[]
            ext=rotated.extent
            oc='%s %s' %(ext.XMin,ext.YMin)
            ya='%s %s' %(ext.XMin,ext.YMax)
            cc='%s %s' %(ext.XMax,ext.YMax)
            arcpy.CreateFishnet_management(fnet, oc,ya, W, L,"","",
                                        "","NO_LABELS", rotated,"POLYGON")
            rects=arcpy.Clip_analysis(fnet, rotated, g)
            for chop in rects:
                if chop.area<A:continue
                FNET.append(chop)
            return FNET
        
        # not sure what's going on here...
        g=arcpy.Geometry()
        PGON=arcpy.CopyFeatures_management(infc,g)[0]
        theList=[PGON];bigList=[]

        nBefore=0
        while True:
            for toCut in theList:
                arcpy.AddMessage("'toCut' variable is {}".format(toCut))
                ## FIND rotation to maximise complete rectangles
                nMax=0
                cPoint=toCut.centroid
                for i in range(36):
                    angle=5*i
                    rotated=ShapeMake(toCut,angle)
                    squares=fnetMake()
                    N=len(squares)
                    if N<=nMax:continue
                    nMax=N
                    keepers=squares[:]
                    bestAngle=angle
                if nMax==0:continue
                arcpy.AddMessage("%s cell(s) found so far" %nMax)
                for item in keepers:
                    rotated=ShapeMake(item,-bestAngle)
                    bigList.append(rotated)
            if nBefore==len(bigList):break
            nBefore=len(bigList)
            arcpy.Erase_analysis(PGON, bigList, erased)
            theList=arcpy.MultipartToSinglepart_management(erased, g)
        arcpy.CopyFeatures_management(bigList,outFC)


        return