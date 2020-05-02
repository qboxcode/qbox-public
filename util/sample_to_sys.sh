#!/bin/bash
# sample_to_sys.sh: extract atomic positions from a sample and generate
# a Qbox input file to create atoms at these positions
# use: sample_to_sys.sh sample.xml > sample.sys
#
xsltproc - $1 << EOF
<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0">
<xsl:output method="text" indent="yes"/>
<xsl:strip-space elements="*"/>
<xsl:template match="/fpmd:sample">
  <xsl:apply-templates/>
</xsl:template>
<xsl:template match="atomset">
  <xsl:apply-templates/>
</xsl:template>
<xsl:template match="unit_cell">
  <xsl:text>set cell </xsl:text>
  <xsl:value-of select="@a"/> <xsl:text> </xsl:text>
  <xsl:value-of select="@b"/> <xsl:text> </xsl:text>
  <xsl:value-of select="@c"/>  <xsl:text>
</xsl:text>
</xsl:template>
<xsl:template match="species">
  <xsl:text>species </xsl:text>
  <xsl:value-of select="@name"/> <xsl:text> </xsl:text>
  <xsl:value-of select="@href"/> <xsl:text>
</xsl:text>
</xsl:template>
<xsl:template match="atom">
  <xsl:text>atom </xsl:text>
  <xsl:value-of select="@name"/> <xsl:text> </xsl:text>
  <xsl:value-of select="@species"/> <xsl:text> </xsl:text>
  <xsl:value-of select="position"/>  <xsl:text> </xsl:text>
  <xsl:value-of select="velocity"/> <xsl:text>
</xsl:text>
</xsl:template>
<xsl:template match="*"/>
</xsl:stylesheet>
EOF
