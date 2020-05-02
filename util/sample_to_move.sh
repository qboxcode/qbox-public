#!/bin/bash
# sample_to_move.sh
# generate a set of Qbox commands to move atoms to the positions
# in a given sample.
# use: sample_to_move.sh sample.xml
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
<xsl:template match="atom">
  <xsl:text>move </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text> to </xsl:text>
  <xsl:value-of select="position"/> <xsl:text>
</xsl:text>
</xsl:template>
<xsl:template match="*"/>
</xsl:stylesheet>
EOF
