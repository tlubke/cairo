/* Copyright © 2022 Matthias Clasen
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA 02110-1335, USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 *
 * The Original Code is the cairo graphics library.
 *
 * Contributor(s):
 *      Matthias Clasen <mclasen@redhat.com>
 */

#include "cairoint.h"
#include "cairo-array-private.h"
#include "cairo-ft-private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if HAVE_FT_GET_COLOR_GLYPH_PAINT

#include <ft2build.h>
#include FT_CONFIG_OPTIONS_H
#include FT_COLOR_H
#include FT_GLYPH_H
#include FT_OUTLINE_H
#include FT_SIZES_H

/* #define DEBUG_COLR 1 */


/* {{{ Utilities */

typedef struct {
  float x, y;
} Point;

static inline float
f16_16 (FT_Fixed f)
{
  return f / (float) (1 << 16);
}

static inline float
f26_6 (FT_F26Dot6 f)
{
  return f / (float) (1 << 6);
}

static inline float
f2_14 (FT_F2Dot14 f)
{
  return f / (float) (1 << 14);
}

static inline float
interpolate (float f0, float f1, float f)
{
  return f0 + f * (f1 - f0);
}

static inline void
interpolate_points (Point *p0, Point *p1, float f, Point *out)
{
  out->x = interpolate (p0->x, p1->x, f);
  out->y = interpolate (p0->y, p1->y, f);
}

static inline void
interpolate_colors (cairo_color_t *c0, cairo_color_t *c1, float f, cairo_color_t *out)
{
  out->red = interpolate (c0->red, c1->red, f);
  out->green = interpolate (c0->green, c1->green, f);
  out->blue = interpolate (c0->blue, c1->blue, f);
  out->alpha = interpolate (c0->alpha, c1->alpha, f);
}

static inline float
dot (Point p, Point q)
{
  return p.x * q.x + p.y * q.y;
}

static inline Point
normalize (Point p)
{
  float len = sqrt (dot (p, p));

  return (Point) { p.x / len, p.y / len };
}

static inline Point
sum (Point p, Point q)
{
  return (Point) { p.x + q.x, p.y + q.y };
}

static inline Point
difference (Point p, Point q)
{
  return (Point) { p.x - q.x, p.y - q.y };
}

static inline Point
scale (Point p, float f)
{
  return (Point) { p.x * f, p.y * f };
}

static cairo_operator_t
cairo_operator (FT_Composite_Mode mode)
{
  switch (mode)
    {
    case FT_COLR_COMPOSITE_CLEAR: return CAIRO_OPERATOR_CLEAR;
    case FT_COLR_COMPOSITE_SRC: return CAIRO_OPERATOR_SOURCE;
    case FT_COLR_COMPOSITE_DEST: return CAIRO_OPERATOR_DEST;
    case FT_COLR_COMPOSITE_SRC_OVER: return CAIRO_OPERATOR_OVER;
    case FT_COLR_COMPOSITE_DEST_OVER: return CAIRO_OPERATOR_DEST_OVER;
    case FT_COLR_COMPOSITE_SRC_IN: return CAIRO_OPERATOR_IN;
    case FT_COLR_COMPOSITE_DEST_IN: return CAIRO_OPERATOR_DEST_IN;
    case FT_COLR_COMPOSITE_SRC_OUT: return CAIRO_OPERATOR_OUT;
    case FT_COLR_COMPOSITE_DEST_OUT: return CAIRO_OPERATOR_DEST_OUT;
    case FT_COLR_COMPOSITE_SRC_ATOP: return CAIRO_OPERATOR_ATOP;
    case FT_COLR_COMPOSITE_DEST_ATOP: return CAIRO_OPERATOR_DEST_ATOP;
    case FT_COLR_COMPOSITE_XOR: return CAIRO_OPERATOR_XOR;
    case FT_COLR_COMPOSITE_PLUS: return CAIRO_OPERATOR_ADD;
    case FT_COLR_COMPOSITE_SCREEN: return CAIRO_OPERATOR_SCREEN;
    case FT_COLR_COMPOSITE_OVERLAY: return CAIRO_OPERATOR_OVERLAY;
    case FT_COLR_COMPOSITE_DARKEN: return CAIRO_OPERATOR_DARKEN;
    case FT_COLR_COMPOSITE_LIGHTEN: return CAIRO_OPERATOR_LIGHTEN;
    case FT_COLR_COMPOSITE_COLOR_DODGE: return CAIRO_OPERATOR_COLOR_DODGE;
    case FT_COLR_COMPOSITE_COLOR_BURN: return CAIRO_OPERATOR_COLOR_BURN;
    case FT_COLR_COMPOSITE_HARD_LIGHT: return CAIRO_OPERATOR_HARD_LIGHT;
    case FT_COLR_COMPOSITE_SOFT_LIGHT: return CAIRO_OPERATOR_SOFT_LIGHT;
    case FT_COLR_COMPOSITE_DIFFERENCE: return CAIRO_OPERATOR_DIFFERENCE;
    case FT_COLR_COMPOSITE_EXCLUSION: return CAIRO_OPERATOR_EXCLUSION;
    case FT_COLR_COMPOSITE_MULTIPLY: return CAIRO_OPERATOR_MULTIPLY;
    case FT_COLR_COMPOSITE_HSL_HUE: return CAIRO_OPERATOR_HSL_HUE;
    case FT_COLR_COMPOSITE_HSL_SATURATION: return CAIRO_OPERATOR_HSL_SATURATION;
    case FT_COLR_COMPOSITE_HSL_COLOR: return CAIRO_OPERATOR_HSL_COLOR;
    case FT_COLR_COMPOSITE_HSL_LUMINOSITY: return CAIRO_OPERATOR_HSL_LUMINOSITY;
    case FT_COLR_COMPOSITE_MAX:
    default:
      assert (0);
    }
}

static cairo_extend_t
cairo_extend (FT_PaintExtend extend)
{
  switch (extend)
    {
    case FT_COLR_PAINT_EXTEND_PAD: return CAIRO_EXTEND_PAD;
    case FT_COLR_PAINT_EXTEND_REPEAT: return CAIRO_EXTEND_REPEAT;
    case FT_COLR_PAINT_EXTEND_REFLECT: return CAIRO_EXTEND_REFLECT;
    default:
      assert (0);
    }
}

/* }}} */
/* {{{ Paths */

typedef struct {
  cairo_array_t points;
  int needs_move;
  int last_move_op;
} path_data_t;

static cairo_status_t
close_path (path_data_t *pdata)
{
  cairo_path_data_t p;
  int status;

  if (pdata->last_move_op < 0)
    return CAIRO_STATUS_SUCCESS;

  p.header.type = CAIRO_PATH_LINE_TO;
  p.header.length = 2;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p = *(cairo_path_data_t *) _cairo_array_index (&pdata->points, pdata->last_move_op + 1);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.header.type = CAIRO_PATH_CLOSE_PATH;
  p.header.length = 1;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  pdata->needs_move = 1;
  pdata->last_move_op = -1;

  return CAIRO_STATUS_SUCCESS;
}

static int
move_to (const FT_Vector *to,
         void *data)
{
  path_data_t *pdata = data;
  cairo_path_data_t p;
  cairo_status_t status;

  status = close_path (pdata);
  if (unlikely (status))
    return status;

  pdata->needs_move = 0;
  pdata->last_move_op = pdata->points.num_elements;

  p.header.type = CAIRO_PATH_MOVE_TO;
  p.header.length = 2;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = to->x / (float) (1 << 6);
  p.point.y = to->y / (float) (1 << 6);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  return CAIRO_STATUS_SUCCESS;
}

static int
line_to (const FT_Vector *to,
         void *data)
{
  path_data_t *pdata = data;
  cairo_path_data_t p;
  cairo_status_t status;

  if (pdata->needs_move)
    return move_to (to, data);

  p.header.type = CAIRO_PATH_LINE_TO;
  p.header.length = 2;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = to->x / (float) (1 << 6);
  p.point.y = to->y / (float) (1 << 6);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  return CAIRO_STATUS_SUCCESS;
}

static int
conic_to (const FT_Vector *control,
          const FT_Vector *to,
          void *data)
{
  path_data_t *pdata = data;
  cairo_path_data_t p;
  double cx, cy;
  double x0, y0;
  double x1, y1;
  double x2, y2;
  double x3, y3;
  cairo_status_t status;

  assert (!pdata->needs_move);
  assert (pdata->points.num_elements > 0);

  p.header.type = CAIRO_PATH_CURVE_TO;
  p.header.length = 4;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p = *(cairo_path_data_t *) _cairo_array_index (&pdata->points, pdata->points.num_elements - 2);

  x0 = p.point.x;
  y0 = p.point.y;

  cx = control->x / (float) (1 << 6);
  cy = control->y / (float) (1 << 6);

  x3 = to->x / (float) (1 << 6);
  y3 = to->y / (float) (1 << 6);

  x1 = x0 + (2./3.) * (cx - x0);
  y1 = y0 + (2./3.) * (cy - y0);

  x2 = x3 + (2./3.) * (cx - x3);
  y2 = y3 + (2./3.) * (cy - y3);

  p.point.x = x1;
  p.point.y = y1;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = x2;
  p.point.y = y2;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = x3;
  p.point.y = y3;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  return CAIRO_STATUS_SUCCESS;
}

static int
cubic_to (const FT_Vector *control1,
          const FT_Vector *control2,
          const FT_Vector *to,
          void *data)
{
  path_data_t *pdata = data;
  cairo_path_data_t p;
  cairo_status_t status;

  assert (!pdata->needs_move);
  assert (pdata->points.num_elements > 0);

  p.header.type = CAIRO_PATH_CURVE_TO;
  p.header.length = 4;
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = control1->x / (float) (1 << 6);
  p.point.y = control1->y / (float) (1 << 6);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = control2->x / (float) (1 << 6);
  p.point.y = control2->y / (float) (1 << 6);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  p.point.x = to->x / (float) (1 << 6);
  p.point.y = to->y / (float) (1 << 6);
  status = _cairo_array_append (&pdata->points, &p);
  if (unlikely (status))
    return status;

  return CAIRO_STATUS_SUCCESS;
}

static cairo_status_t
get_path_for_glyph (FT_Face face,
                    unsigned int glyph,
                    cairo_path_t **path)
{
  path_data_t pdata;
  FT_Outline_Funcs callbacks = {
    move_to, line_to, conic_to, cubic_to, 0, 0
  };
  FT_Error error;

  *path = NULL;

  error = FT_Load_Glyph (face, glyph, FT_LOAD_DEFAULT);
  if (error != 0)
    return CAIRO_STATUS_INVALID_PATH_DATA;

  _cairo_array_init (&pdata.points, sizeof (cairo_path_data_t));
  pdata.needs_move = 1;
  pdata.last_move_op = -1;

  if (FT_Outline_Decompose (&face->glyph->outline, &callbacks, &pdata) != 0)
    {
      _cairo_array_fini (&pdata.points);
      return CAIRO_STATUS_INVALID_PATH_DATA;
    }

  close_path (&pdata);

  *path = _cairo_malloc (sizeof (cairo_path_t));
  if (unlikely (path == NULL))
    return CAIRO_STATUS_NO_MEMORY;

  (*path)->status = CAIRO_STATUS_SUCCESS;
  (*path)->data = (cairo_path_data_t *) pdata.points.elements;
  (*path)->num_data = pdata.points.num_elements;

  return CAIRO_STATUS_SUCCESS;
}


/* }}} */
/* {{{ cairo_colr_glyph_render_t implementation */

typedef struct cairo_colr_glyph_render_t cairo_colr_glyph_render_t;

struct cairo_colr_glyph_render_t {
  FT_Face face;
  cairo_pattern_t *foreground_color;
  FT_Color *palette;
  unsigned int num_palette_entries;
  int level;
};

static cairo_status_t draw_paint (cairo_colr_glyph_render_t *render,
                                  FT_OpaquePaint *paint,
                                  cairo_t *cr);

static cairo_status_t
draw_paint_colr_layers (cairo_colr_glyph_render_t *render,
                        FT_PaintColrLayers *colr_layers,
                        cairo_t *cr)
{
  FT_OpaquePaint paint;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintColrLayers\n", 2 * render->level, "");
#endif

  while (FT_Get_Paint_Layers (render->face, &colr_layers->layer_iterator, &paint))
    {
      cairo_push_group (cr);
      status = draw_paint (render, &paint, cr);
      cairo_pop_group_to_source (cr);
      cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
      cairo_paint (cr);

      if (unlikely (status))
        break;
    }

  return status;
}

static void
get_palette_color (cairo_colr_glyph_render_t *render,
                   FT_ColorIndex *ci,
                   cairo_color_t *color,
                   cairo_bool_t *is_foreground_color)
{
    cairo_bool_t foreground = FALSE;

    if (ci->palette_index == 0xffff || ci->palette_index >= render->num_palette_entries) {
        color->red = 0;
        color->green = 0;
        color->blue = 0;
        foreground = TRUE;
    } else {
        FT_Color c = render->palette[ci->palette_index];
        color->red = c.red / 255.0;
        color->green = c.green / 255.0;
        color->blue = c.blue / 255.0;
    }
    color->alpha = f2_14 (ci->alpha);
    if (foreground)
        *is_foreground_color = TRUE;
}

static cairo_status_t
draw_paint_solid (cairo_colr_glyph_render_t *render,
                  FT_PaintSolid *solid,
                  cairo_t *cr)
{
  cairo_color_t color;
  cairo_bool_t is_foreground_color;

#if DEBUG_COLR
  printf ("%*sDraw PaintSolid\n", 2 * render->level, "");
#endif

  get_palette_color (render, &solid->color, &color, &is_foreground_color);
  if (is_foreground_color) {
      cairo_set_source (cr, render->foreground_color);
  } else
      cairo_set_source_rgba (cr, color.red, color.green, color.blue, color.alpha);

  cairo_paint (cr);

  return CAIRO_STATUS_SUCCESS;
}

typedef struct {
  cairo_color_t color;
  float position;
} ColorStop;

typedef struct {
  int n_stops;
  ColorStop *stops;
} ColorLine;

static void
free_colorline (ColorLine *cl)
{
  free (cl->stops);
  free (cl);
}

static int
compare_stops (const void *p1, const void *p2)
{
  const ColorStop *c1 = p1;
  const ColorStop *c2 = p2;

  if (c1->position < c2->position)
    return -1;
  else if (c1->position > c2->position)
    return 1;
  else
    return 0;
}

static ColorLine *
read_colorline (cairo_colr_glyph_render_t *render,
                FT_ColorLine *colorline)
{
  ColorLine *cl;
  FT_ColorStop stop;
  int i;

  cl = calloc (1, sizeof (ColorLine));
  if (unlikely (cl == NULL))
    return NULL;

  cl->n_stops = colorline->color_stop_iterator.num_color_stops;
  cl->stops = calloc (cl->n_stops, sizeof (ColorStop));
  if (unlikely (cl->stops == NULL)) {
    free (cl);
    return NULL;
  }

  i = 0;
  while (FT_Get_Colorline_Stops (render->face, &stop, &colorline->color_stop_iterator))
    {
      cl->stops[i].position = f16_16 (stop.stop_offset);
      get_palette_color (render, &stop.color, &cl->stops[i].color, NULL);
      i++;
    }

  qsort (cl->stops, cl->n_stops, sizeof (ColorStop), compare_stops);

  return cl;
}

static void
reduce_anchors (FT_PaintLinearGradient *gradient,
                Point *pp0,
                Point *pp1)
{
  Point p0, p1, p2;
  Point q1, q2;
  float s;
  float k;

  p0.x = f16_16 (gradient->p0.x);
  p0.y = f16_16 (gradient->p0.y);
  p1.x = f16_16 (gradient->p1.x);
  p1.y = f16_16 (gradient->p1.y);
  p2.x = f16_16 (gradient->p2.x);
  p2.y = f16_16 (gradient->p2.y);

  q2.x = p2.x - p0.x;
  q2.y = p2.y - p0.y;
  q1.x = p1.x - p0.x;
  q1.y = p1.y - p0.y;

  s = q2.x * q2.x + q2.y * q2.y;
  if (s < 0.000001)
    {
      pp0->x = p0.x; pp0->y = p0.y;
      pp1->x = p1.x; pp1->y = p1.y;
      return;
    }

  k = (q2.x * q1.x + q2.y * q1.y) / s;
  pp0->x = p0.x;
  pp0->y = p0.y;
  pp1->x = p1.x - k * q2.x;
  pp1->y = p1.y - k * q2.y;
}

static void
normalize_colorline (ColorLine *cl,
                     float *out_min,
                     float *out_max)
{
  float min, max;

  *out_min = 0.;
  *out_max = 1.;

  min = max = cl->stops[0].position;
  for (int i = 0; i < cl->n_stops; i++)
    {
      ColorStop *stop = &cl->stops[i];
      min = MIN (min, stop->position);
      max = MAX (max, stop->position);
    }

  if (min != max)
    {
      for (int i = 0; i < cl->n_stops; i++)
        {
          ColorStop *stop = &cl->stops[i];
          stop->position = (stop->position - min) / (max - min);
        }
      *out_min = min;
      *out_max = max;
    }
}

static cairo_status_t
draw_paint_linear_gradient (cairo_colr_glyph_render_t *render,
                            FT_PaintLinearGradient *gradient,
                            cairo_t *cr)
{
  ColorLine *cl;
  Point p0, p1;
  Point pp0, pp1;
  cairo_pattern_t *pattern;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;
  float min, max;

#if DEBUG_COLR
  printf ("%*sDraw PaintLinearGradient\n", 2 * render->level, "");
#endif

  cl = read_colorline (render, &gradient->colorline);
  if (unlikely (cl == NULL))
    return CAIRO_STATUS_NO_MEMORY;

  /* cairo only allows stop positions between 0 and 1 */
  normalize_colorline (cl, &min, &max);
  reduce_anchors (gradient, &p0, &p1);
  interpolate_points (&p0, &p1, min, &pp0);
  interpolate_points (&p0, &p1, max, &pp1);

  pattern = cairo_pattern_create_linear (pp0.x, pp0.y, pp1.x, pp1.y);

  cairo_pattern_set_extend (pattern, cairo_extend (gradient->colorline.extend));

  for (int i = 0; i < cl->n_stops; i++)
    {
      ColorStop *stop = &cl->stops[i];
      cairo_pattern_add_color_stop_rgba (pattern, stop->position,
                                         stop->color.red, stop->color.green, stop->color.blue, stop->color.alpha);
    }

  cairo_set_source (cr, pattern);
  cairo_paint (cr);

  cairo_pattern_destroy (pattern);

  free_colorline (cl);

  return status;
}

static cairo_status_t
draw_paint_radial_gradient (cairo_colr_glyph_render_t *render,
                            FT_PaintRadialGradient *gradient,
                            cairo_t *cr)
{
  ColorLine *cl;
  Point start, end;
  Point start1, end1;
  float start_radius, end_radius;
  float start_radius1, end_radius1;
  float min, max;
  cairo_pattern_t *pattern;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintRadialGradient\n", 2 * render->level, "");
#endif

  cl = read_colorline (render, &gradient->colorline);
  if (unlikely (cl == NULL))
    return CAIRO_STATUS_NO_MEMORY;

  start.x = f16_16 (gradient->c0.x);
  start.y = f16_16 (gradient->c0.y);
  end.x = f16_16 (gradient->c1.x);
  end.y = f16_16 (gradient->c1.y);

  start_radius = f16_16 (gradient->r0);
  end_radius = f16_16 (gradient->r1);

  /* cairo only allows stop positions between 0 and 1 */
  normalize_colorline (cl, &min, &max);
  interpolate_points (&start, &end, min, &start1);
  interpolate_points (&start, &end, max, &end1);
  start_radius1 = interpolate (start_radius, end_radius, min);
  end_radius1 = interpolate (start_radius, end_radius, max);

  pattern = cairo_pattern_create_radial (start1.x, start1.y, start_radius1,
                                         end1.x, end1.y, end_radius1);

  cairo_pattern_set_extend (pattern, cairo_extend (gradient->colorline.extend));

  for (int i = 0; i < cl->n_stops; i++)
    {
      ColorStop *stop = &cl->stops[i];
      cairo_pattern_add_color_stop_rgba (pattern, stop->position,
                                         stop->color.red, stop->color.green, stop->color.blue, stop->color.alpha);
    }

  cairo_set_source (cr, pattern);
  cairo_paint (cr);

  cairo_pattern_destroy (pattern);

  free_colorline (cl);

  return status;
}

typedef struct {
  Point center, p0, c0, c1, p1;
  cairo_color_t color0, color1;
} Patch;

static void
add_patch (cairo_pattern_t *pattern, Point *center, Patch *p)
{
  cairo_mesh_pattern_begin_patch (pattern);
  cairo_mesh_pattern_move_to (pattern, center->x, center->y);
  cairo_mesh_pattern_line_to (pattern, p->p0.x, p->p0.y);
  cairo_mesh_pattern_curve_to (pattern,
                               p->c0.x, p->c0.y,
                               p->c1.x, p->c1.y,
                               p->p1.x, p->p1.y);
  cairo_mesh_pattern_line_to (pattern, center->x, center->y);
  cairo_mesh_pattern_set_corner_color_rgba (pattern, 0,
                                            p->color0.red,
                                            p->color0.green,
                                            p->color0.blue,
                                            p->color0.alpha);
  cairo_mesh_pattern_set_corner_color_rgba (pattern, 1,
                                            p->color0.red,
                                            p->color0.green,
                                            p->color0.blue,
                                            p->color0.alpha);
  cairo_mesh_pattern_set_corner_color_rgba (pattern, 2,
                                            p->color1.red,
                                            p->color1.green,
                                            p->color1.blue,
                                            p->color1.alpha);
  cairo_mesh_pattern_set_corner_color_rgba (pattern, 3,
                                            p->color1.red,
                                            p->color1.green,
                                            p->color1.blue,
                                            p->color1.alpha);
  cairo_mesh_pattern_end_patch (pattern);
}

#define MAX_ANGLE (M_PI / 8.)

static void
add_sweep_gradient_patches1 (Point *center, float radius,
                             float a0, cairo_color_t *c0,
                             float a1, cairo_color_t *c1,
                             cairo_pattern_t *pattern)
{

  int num_splits;
  Point p0;
  cairo_color_t color0, color1;

  num_splits = ceilf (fabs (a1 - a0) / MAX_ANGLE);
  p0 = (Point) { cosf (a0), sinf (a0) };
  color0 = *c0;

  for (int a = 0; a < num_splits; a++)
    {
      float k = (a + 1.) / num_splits;
      float angle1;
      Point p1;
      Point A, U;
      Point C0, C1;
      Patch patch;

      angle1 = interpolate (a0, a1, k);
      interpolate_colors (c0, c1, k, &color1);

      patch.color0 = color0;
      patch.color1 = color1;

      p1 = (Point) { cosf (angle1), sinf (angle1) };
      patch.p0 = sum (*center, scale (p0, radius));
      patch.p1 = sum (*center, scale (p1, radius));

      A = normalize (sum (p0, p1));
      U = (Point) { -A.y, A.x };
      C0 = sum (A, scale (U, dot (difference (p0, A), p0) / dot (U, p0)));
      C1 = sum (A, scale (U, dot (difference (p1, A), p1) / dot (U, p1)));
      patch.c0 = sum (*center, scale (sum (C0, scale (difference (C0, p0), 0.33333)), radius));
      patch.c1 = sum (*center, scale (sum (C1, scale (difference (C1, p1), 0.33333)), radius));

      add_patch (pattern, center, &patch);

      p0 = p1;
      color0 = color1;
    }
}

static void
add_sweep_gradient_patches (ColorLine *cl,
                            cairo_extend_t extend,
                            Point *center,
                            float radius,
                            float start_angle,
                            float end_angle,
                            cairo_pattern_t *pattern)
{
  float *angles;
  cairo_color_t color0, color1;

  if (start_angle == end_angle)
    {
      if (extend == CAIRO_EXTEND_PAD)
        {
          if (start_angle > 0)
            add_sweep_gradient_patches1 (center, radius,
                                         0.,          &cl->stops[0].color,
                                         start_angle, &cl->stops[0].color,
                                         pattern);
          if (end_angle < 2 * M_PI)
            add_sweep_gradient_patches1 (center, radius,
                                         end_angle, &cl->stops[cl->n_stops - 1].color,
                                         2 * M_PI,  &cl->stops[cl->n_stops - 1].color,
                                         pattern);
        }
      return;
    }

  assert (start_angle != end_angle);

  angles = alloca (sizeof (float) * cl->n_stops);

  for (int i = 0; i < cl->n_stops; i++)
    angles[i] = start_angle + cl->stops[i].position * (end_angle - start_angle);

  /* handle directions */
  if (end_angle < start_angle)
    {
      for (int i = 0; i < cl->n_stops - 1 - i; i++)
        {
          ColorStop stop = cl->stops[i];
          float a = angles[i];
          cl->stops[i] = cl->stops[cl->n_stops - 1 - i];
          cl->stops[cl->n_stops - 1 - i] = stop;
          angles[i] = angles[cl->n_stops - 1 - i];
          angles[cl->n_stops - 1 - i] = a;
        }
    }

  if (extend == CAIRO_EXTEND_PAD)
    {
      int pos;

      color0 = cl->stops[0].color;
      for (pos = 0; pos < cl->n_stops; pos++)
        {
          if (angles[pos] >= 0)
            {
              if (pos > 0)
                {
                  float k = (0 - angles[pos - 1]) / (angles[pos] - angles[pos - 1]);
                  interpolate_colors (&cl->stops[pos - 1].color, &cl->stops[pos].color, k, &color0);
                }
              break;
            }
        }
      if (pos == cl->n_stops)
        {
          /* everything is below 0 */
          color0 = cl->stops[cl->n_stops - 1].color;
          add_sweep_gradient_patches1 (center, radius,
                                       0.,       &color0,
                                       2 * M_PI, &color0,
                                       pattern);
          return;
        }

      add_sweep_gradient_patches1 (center, radius,
                                   0.,          &color0,
                                   angles[pos], &cl->stops[pos].color,
                                   pattern);

      for (pos++; pos < cl->n_stops; pos++)
        {
          if (angles[pos] <= 2 * M_PI)
            {
              add_sweep_gradient_patches1 (center, radius,
                                           angles[pos - 1], &cl->stops[pos - 1].color,
                                           angles[pos],     &cl->stops[pos].color,
                                           pattern);
            }
          else
            {
              float k = (2 * M_PI - angles[pos - 1]) / (angles[pos] - angles[pos - 1]);
              interpolate_colors (&cl->stops[pos - 1].color, &cl->stops[pos].color, k, &color1);
              add_sweep_gradient_patches1 (center, radius,
                                           angles[pos - 1], &cl->stops[pos - 1].color,
                                           2 * M_PI,        &color1,
                                           pattern);
              break;
            }
        }

      if (pos == cl->n_stops)
        {
          /* everything is below 2*M_PI */
          color0 = cl->stops[cl->n_stops - 1].color;
          add_sweep_gradient_patches1 (center, radius,
                                       angles[cl->n_stops - 1], &color0,
                                       2 * M_PI,                &color0,
                                       pattern);
          return;
        }
    }
  else
    {
      int k;
      float span;

      span = angles[cl->n_stops - 1] - angles[0];
      k = 0;
      if (angles[0] >= 0)
        {
          float ss = angles[0];
          while (ss > 0)
            {
              if (span > 0)
                {
                  ss -= span;
                  k--;
                }
              else
                {
                  ss += span;
                  k++;
                }
            }
        }
      else if (angles[0] < 0)
        {
          float ee = angles[cl->n_stops - 1];
          while (ee < 0)
            {
              if (span > 0)
                {
                  ee += span;
                  k++;
                }
              else
                {
                  ee -= span;
                  k--;
                }
            }
        }

      //assert (angles[0] + k * span <= 0 && 0 < angles[cl->n_stops - 1] + k * span);

      for (int l = k; TRUE; l++)
        {
          for (int i = 1; i < cl->n_stops; i++)
            {
              float a0, a1;
              cairo_color_t *c0, *c1;

              if ((l % 2 != 0) && (extend == CAIRO_EXTEND_REFLECT))
                {
                  a0 = angles[0] + angles[cl->n_stops - 1] - angles[cl->n_stops - 1 - (i-1)] + l * span;
                  a1 = angles[0] + angles[cl->n_stops - 1] - angles[cl->n_stops - 1 - i] + l * span;
                  c0 = &cl->stops[cl->n_stops - 1 - (i-1)].color;
                  c1 = &cl->stops[cl->n_stops - 1 - i].color;
                }
              else
                {
                  a0 = angles[i-1] + l * span;
                  a1 = angles[i] + l * span;
                  c0 = &cl->stops[i-1].color;
                  c1 = &cl->stops[i].color;
                }

              if (a1 < 0)
                continue;
              if (a0 < 0)
                {
                  cairo_color_t color;
                  float f = (0 - a0)/(a1 - a0);
                  interpolate_colors (c0, c1, f, &color);
                  add_sweep_gradient_patches1 (center, radius,
                                               0,  &color,
                                               a1, c1,
                                               pattern);
                }
              else if (a1 >= 2 * M_PI)
                {
                  cairo_color_t color;
                  float f = (2 * M_PI - a0)/(a1 - a0);
                  interpolate_colors (c0, c1, f, &color);
                  add_sweep_gradient_patches1 (center, radius,
                                               a0,       c0,
                                               2 * M_PI, &color,
                                               pattern);
                  return;
                }
              else
                {
                  add_sweep_gradient_patches1 (center, radius,
                                               a0, c0,
                                               a1, c1,
                                               pattern);
                }
            }
        }
    }
}

static cairo_status_t
draw_paint_sweep_gradient (cairo_colr_glyph_render_t *render,
                           FT_PaintSweepGradient *gradient,
                           cairo_t *cr)
{
  ColorLine *cl;
  Point center;
  float start_angle, end_angle;
  double x1, y1, x2, y2;
  double max_x, max_y, R;
  cairo_pattern_t *pattern;
  cairo_extend_t extend;

#if DEBUG_COLR
  printf ("%*sDraw PaintSweepGradient\n", 2 * render->level, "");
#endif

  cl = read_colorline (render, &gradient->colorline);
  if (unlikely (cl == NULL))
    return CAIRO_STATUS_NO_MEMORY;

  center.x = f16_16 (gradient->center.x);
  center.y = f16_16 (gradient->center.y);
  start_angle = (f16_16 (gradient->start_angle) + 1) * M_PI;
  end_angle = (f16_16 (gradient->end_angle) + 1) * M_PI;

  pattern = cairo_pattern_create_mesh ();

  cairo_clip_extents (cr, &x1, &y1, &x2, &y2);
  max_x = MAX ((x1 - center.x) * (x1 - center.x), (x2 - center.x) * (x2 - center.x));
  max_y = MAX ((y1 - center.y) * (y1 - center.y), (y2 - center.y) * (y2 - center.y));
  R = sqrt (max_x + max_y);

  extend = cairo_extend (gradient->colorline.extend);

  add_sweep_gradient_patches (cl, extend, &center, R, start_angle, end_angle, pattern);

  cairo_set_source (cr, pattern);
  cairo_paint (cr);

  cairo_pattern_destroy (pattern);

  free_colorline (cl);

  return CAIRO_STATUS_SUCCESS;
}

static cairo_status_t
draw_paint_glyph (cairo_colr_glyph_render_t *render,
                  FT_PaintGlyph *glyph,
                  cairo_t *cr)
{
  cairo_path_t *path;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintGlyph\n", 2 * render->level, "");
#endif

  status = get_path_for_glyph (render->face, glyph->glyphID, &path);
  if (unlikely (status))
    goto cleanup;

  cairo_save (cr);

  cairo_new_path (cr);
  cairo_append_path (cr, path);
  cairo_clip (cr);

  status = draw_paint (render, &glyph->paint, cr);

  cairo_restore (cr);

cleanup:

  if (path)
    cairo_path_destroy (path);

  return status;
}

static cairo_status_t draw_colr_glyph (cairo_colr_glyph_render_t *render,
                                       unsigned int glyph,
                                       FT_Color_Root_Transform  root,
                                       cairo_t *cr);

static cairo_status_t
draw_paint_colr_glyph (cairo_colr_glyph_render_t *render,
                       FT_PaintColrGlyph *colr_glyph,
                       cairo_t *cr)
{
  //printf ("%*sDraw PaintColrGlyph\n", 2 * render->level, "");

  return draw_colr_glyph (render, colr_glyph->glyphID, FT_COLOR_NO_ROOT_TRANSFORM, cr);
}

static cairo_status_t
draw_paint_transform (cairo_colr_glyph_render_t *render,
                      FT_PaintTransform *transform,
                      cairo_t *cr)
{
  cairo_matrix_t t;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintTransform\n", 2 * render->level, "");
#endif

  cairo_matrix_init (&t,
                     f16_16 (transform->affine.xx),
                     f16_16 (transform->affine.yx),
                     f16_16 (transform->affine.xy),
                     f16_16 (transform->affine.yy),
                     f16_16 (transform->affine.dx),
                     f16_16 (transform->affine.dy));

  cairo_save (cr);

  cairo_transform (cr, &t);
  status = draw_paint (render, &transform->paint, cr);

  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint_translate (cairo_colr_glyph_render_t *render,
                      FT_PaintTranslate *translate,
                      cairo_t *cr)
{
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintTranslate\n", 2 * render->level, "");
#endif

  cairo_save (cr);

  cairo_translate (cr, f16_16 (translate->dx), f16_16 (translate->dy));
  status = draw_paint (render, &translate->paint, cr);

  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint_rotate (cairo_colr_glyph_render_t *render,
                   FT_PaintRotate *rotate,
                   cairo_t *cr)
{
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintRotate\n", 2 * render->level, "");
#endif

  cairo_save (cr);

  cairo_translate (cr, f16_16 (rotate->center_x), f16_16 (rotate->center_y));
  cairo_rotate (cr, f16_16 (rotate->angle) * M_PI);
  cairo_translate (cr, - f16_16 (rotate->center_x), - f16_16 (rotate->center_y));
  status = draw_paint (render, &rotate->paint, cr);

  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint_scale (cairo_colr_glyph_render_t *render,
                  FT_PaintScale *scale,
                  cairo_t *cr)
{
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintScale\n", 2 * render->level, "");
#endif

  cairo_save (cr);

  cairo_translate (cr, f16_16 (scale->center_x), f16_16 (scale->center_y));
  cairo_scale (cr, f16_16 (scale->scale_x), f16_16 (scale->scale_y));
  cairo_translate (cr, - f16_16 (scale->center_x), - f16_16 (scale->center_y));
  status = draw_paint (render, &scale->paint, cr);

  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint_skew (cairo_colr_glyph_render_t *render,
                 FT_PaintSkew *skew,
                 cairo_t *cr)
{
  cairo_matrix_t s;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

#if DEBUG_COLR
  printf ("%*sDraw PaintSkew\n", 2 * render->level, "");
#endif

  cairo_save (cr);

  cairo_translate (cr, f16_16 (skew->center_x), f16_16 (skew->center_y));
  cairo_matrix_init (&s, 1., tan (f16_16 (skew->y_skew_angle) * M_PI), - tan (f16_16 (skew->x_skew_angle) * M_PI), 1., 0., 0.);
  cairo_transform (cr, &s);
  cairo_translate (cr, - f16_16 (skew->center_x), - f16_16 (skew->center_y));
  status = draw_paint (render, &skew->paint, cr);

  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint_composite (cairo_colr_glyph_render_t *render,
                      FT_PaintComposite *composite,
                      cairo_t *cr)
{
  cairo_status_t status = CAIRO_STATUS_SUCCESS;


#if DEBUG_COLR
  printf ("%*sDraw PaintComposite\n", 2 * render->level, "");
#endif

  cairo_save (cr);

  cairo_push_group (cr);
  status = draw_paint (render, &composite->backdrop_paint, cr);
  if (unlikely (status)) {
    cairo_pattern_destroy (cairo_pop_group (cr));
    goto cleanup;
  }

  cairo_push_group (cr);
  status = draw_paint (render, &composite->source_paint, cr);
  if (unlikely (status)) {
    cairo_pattern_destroy (cairo_pop_group (cr));
    cairo_pattern_destroy (cairo_pop_group (cr));
    goto cleanup;
  }

  cairo_pop_group_to_source (cr);
  cairo_set_operator (cr, cairo_operator (composite->composite_mode));
  cairo_paint (cr);
  cairo_pop_group_to_source (cr);
  cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
  cairo_paint (cr);

cleanup:
  cairo_restore (cr);

  return status;
}

static cairo_status_t
draw_paint (cairo_colr_glyph_render_t *render,
            FT_OpaquePaint *paint,
            cairo_t *cr)
{
  FT_COLR_Paint p;
  FT_Size orig_size;
  FT_Size unscaled_size;
  FT_Matrix orig_transform;
  FT_Vector orig_delta;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

  assert (cairo_status (cr) == CAIRO_STATUS_SUCCESS);

  if (!FT_Get_Paint (render->face, *paint, &p))
    return CAIRO_STATUS_NO_MEMORY;

  if (render->level == 0)
    {
      /* Now that the FT_Get_Paint call has applied the root transform,
       * make the face unscaled and untransformed, so we can load glyph
       * contours.
       */

      FT_Matrix transform;
      FT_Vector delta;

      orig_size = render->face->size;
      FT_New_Size (render->face, &unscaled_size);
      FT_Activate_Size (unscaled_size);
      FT_Set_Char_Size (render->face, render->face->units_per_EM << 6, 0, 0, 0);

      transform.xx = transform.yy = 1 << 16;
      transform.xy = transform.yx = 0;
      delta.x = delta.y = 0;

      FT_Get_Transform (render->face, &orig_transform, &orig_delta);
      FT_Set_Transform (render->face, &transform, &delta);
    }

  render->level++;

  switch (p.format)
    {
    case FT_COLR_PAINTFORMAT_COLR_LAYERS:
      status = draw_paint_colr_layers (render, &p.u.colr_layers, cr);
      break;
    case FT_COLR_PAINTFORMAT_SOLID:
      status = draw_paint_solid (render, &p.u.solid, cr);
      break;
    case FT_COLR_PAINTFORMAT_LINEAR_GRADIENT:
      status = draw_paint_linear_gradient (render, &p.u.linear_gradient, cr);
      break;
    case FT_COLR_PAINTFORMAT_RADIAL_GRADIENT:
      status = draw_paint_radial_gradient (render, &p.u.radial_gradient, cr);
      break;
    case FT_COLR_PAINTFORMAT_SWEEP_GRADIENT:
      status = draw_paint_sweep_gradient (render, &p.u.sweep_gradient, cr);
      break;
    case FT_COLR_PAINTFORMAT_GLYPH:
      status = draw_paint_glyph (render, &p.u.glyph, cr);
      break;
    case FT_COLR_PAINTFORMAT_COLR_GLYPH:
      status = draw_paint_colr_glyph (render, &p.u.colr_glyph, cr);
      break;
    case FT_COLR_PAINTFORMAT_TRANSFORM:
      status = draw_paint_transform (render, &p.u.transform, cr);
      break;
    case FT_COLR_PAINTFORMAT_TRANSLATE:
      status = draw_paint_translate (render, &p.u.translate, cr);
      break;
    case FT_COLR_PAINTFORMAT_ROTATE:
      status = draw_paint_rotate (render, &p.u.rotate, cr);
      break;
    case FT_COLR_PAINTFORMAT_SCALE:
      status = draw_paint_scale (render, &p.u.scale, cr);
      break;
    case FT_COLR_PAINTFORMAT_SKEW:
      status = draw_paint_skew (render, &p.u.skew, cr);
      break;
    case FT_COLR_PAINTFORMAT_COMPOSITE:
      status = draw_paint_composite (render, &p.u.composite, cr);
      break;
    case FT_COLR_PAINT_FORMAT_MAX:
    case FT_COLR_PAINTFORMAT_UNSUPPORTED:
    default:
      assert (0);
    }

  render->level--;

  if (render->level == 0)
    {
      FT_Set_Transform (render->face, &orig_transform, &orig_delta);
      FT_Activate_Size (orig_size);
      FT_Done_Size (unscaled_size);
    }

  return status;
}

static cairo_status_t
draw_colr_glyph (cairo_colr_glyph_render_t *render,
                 unsigned int glyph,
                 FT_Color_Root_Transform  root,
                 cairo_t *cr)
{
  FT_OpaquePaint paint = { NULL, 0 };
  FT_ClipBox box;
  cairo_status_t status = CAIRO_STATUS_SUCCESS;

  cairo_save (cr);

  if (FT_Get_Color_Glyph_ClipBox (render->face, glyph, &box))
    {
      float xmin, ymin, xmax, ymax;

      xmin = f26_6 (box.bottom_left.x);
      ymin = f26_6 (box.bottom_left.y);
      xmax = f26_6 (box.top_right.x);
      ymax = f26_6 (box.top_right.y);

      cairo_new_path (cr);
      cairo_rectangle (cr, xmin, ymin, xmax - xmin, ymax - ymin);
      cairo_clip (cr);
      printf("%f %f %f %f\n", xmin, ymin, xmax, ymax);
    }

  if (FT_Get_Color_Glyph_Paint (render->face, glyph, root, &paint))
    {
      status = draw_paint (render, &paint, cr);
    }

  cairo_restore (cr);

  return status;
}

/* }}} */
/* {{{ cairo_colr_glyph_render_t API */

/* Create an image surface and render the glyph onto it,
 * using the given colors.
 */
cairo_status_t
_cairo_render_colr_v1_glyph (FT_Face face,
                             unsigned long glyph,
                             FT_UShort palette_index,
                             cairo_t              *cr)
{
  cairo_status_t status = CAIRO_STATUS_SUCCESS;
  cairo_colr_glyph_render_t *colr_render = NULL;
  FT_Color *palette = NULL;
  FT_Palette_Data palette_data;

#if DEBUG_COLR
  printf ("_cairo_render_colr_glyph  glyph index: %ld\n", glyph);
#endif

  colr_render = _cairo_malloc (sizeof (cairo_colr_glyph_render_t));
  if (unlikely (colr_render == NULL)) {
    status = _cairo_error (CAIRO_STATUS_NO_MEMORY);
    goto cleanup;
  }

  if (FT_Palette_Data_Get (face, &palette_data) == 0 && palette_data.num_palettes > 0) {
    if (palette_index >= palette_data.num_palettes)
      palette_index = CAIRO_COLOR_PALETTE_DEFAULT;
    if (FT_Palette_Select (face, palette_index, &palette) != 0)
      palette = NULL;
  }

  colr_render->face = face;
  colr_render->palette = palette;
  colr_render->num_palette_entries = palette_data.num_palette_entries;
  colr_render->foreground_color = cairo_pattern_reference (cairo_get_source (cr));
  colr_render->level = 0;

  status = draw_colr_glyph (colr_render,
                            glyph,
                            FT_COLOR_INCLUDE_ROOT_TRANSFORM,
                            cr);
  
cleanup:

  cairo_pattern_destroy (colr_render->foreground_color);

  if (colr_render)
    free (colr_render);

  return status;
}

/* }}} */

#endif /* HAVE_FT_GET_COLOR_GLYPH_PAINT */

/* vim:set foldmethod=marker expandtab: */
