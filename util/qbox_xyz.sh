#!/bin/bash
# qbox_xyz.sh: get atomic positions from an MD simulation in xyz format
# use: qbox_xyz.sh  mdrun.r > file.xyz
#
xsltproc - $1 << EOF
<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
xmlns:fpmd="http://www.quantum-simulation.org/ns/fpmd/fpmd-1.0">
<xsl:output method="text" indent="yes"/>
<xsl:strip-space elements="*"/>
<xsl:template match="/">
  <xsl:apply-templates select="//iteration/atomset"/>
</xsl:template>

<xsl:template match="iteration/atomset">
  <xsl:value-of select="count(child::atom)"/> <xsl:text>
</xsl:text>
<xsl:number count="iteration"/> <xsl:text>
</xsl:text>
  <xsl:apply-templates select="atom"/>
</xsl:template>

<xsl:template match="atom">
  <xsl:variable name="sym" select="substring(@name,1,2)"/>
  <xsl:variable name="symbol" select="translate(\$sym,'0123456789_-:.',' ')"/>
  <xsl:value-of select="\$symbol"/> <xsl:text> </xsl:text>
  <xsl:variable name="pos" select="normalize-space(position)"/>
  <xsl:variable name="x" select="substring-before(\$pos,' ')"/>
  <xsl:variable name="y" select="substring-before(substring-after(\$pos,' '),' ')"/>
  <xsl:variable name="z" select="substring-after(substring-after(\$pos,' '),' ')"/>
  <xsl:value-of select="\$x * 0.529177"/> <xsl:text> </xsl:text>
  <xsl:value-of select="\$y * 0.529177"/> <xsl:text> </xsl:text>
  <xsl:value-of select="\$z * 0.529177"/> <xsl:text>
</xsl:text>
</xsl:template>

<xsl:template match="*"/>
</xsl:stylesheet>
EOF
