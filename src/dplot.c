#include "dplot.h"

void kh_ins(kh_t *h, uint64_t n, uint64_t v)
{
	int absent;
	khint_t k;
	if ((k = kh_get(h, n)) == kh_end(h))
	{
		k = kh_put(h, n, &absent);
		kh_val(h, k) = v;
	}
	else
		kh_val(h, k) += v;
}

uint64_t kh_xval(const kh_t *h, const uint64_t n)
{
	khint_t k;
	k = kh_get(h, n);
	if (k == kh_end(h))
		return 0;
	else
		return kh_val(h, k);
}

void swap(char *xp, char *yp)
{
	*xp = *xp ^ *yp;
	*yp = *xp ^ *yp;
	*xp = *xp ^ *yp;
}
void reverse(char *str)
{
	int i = 0, j = strlen(str) - 1;
	while (i < j)
		swap(str + i++, str + j--);
}

void draw_rrect(cairo_t *cr)
{
	// a custom shape that could be wrapped in a function
	double x         = 0,        // parameters like cairo_rectangle
	       y         = 0,
	       width         = WIDTH,
	       height        = HEIGHT,
	       aspect        = 1.0,     // aspect ratio
	       corner_radius = height / 60.0;   // and corner curvature radius
	double radius = corner_radius / aspect;
	double degrees = M_PI / 180.0;
	cairo_new_sub_path (cr);
	cairo_arc (cr, x + width - radius, y + radius, radius, -90 * degrees, 0 * degrees);
	cairo_arc (cr, x + width - radius, y + height - radius, radius, 0 * degrees, 90 * degrees);
	cairo_arc (cr, x + radius, y + height - radius, radius, 90 * degrees, 180 * degrees);
	cairo_arc (cr, x + radius, y + radius, radius, 180 * degrees, 270 * degrees);
	cairo_close_path (cr);
	cairo_set_source_rgba (cr, .96, .96, .96, .36);
	cairo_fill(cr);
}

void draw_arrow(cairo_t *cr, double start_x, double start_y, double end_x, double end_y)
{
	double angle = atan2(end_y - start_y, end_x - start_x) + M_PI;
	double arrow_degrees_ = M_PI / 15;
	double arrow_length_ = 15;
	double x1 = end_x + arrow_length_ * cos(angle - arrow_degrees_);
	double y1 = end_y + arrow_length_ * sin(angle - arrow_degrees_);
	double x2 = end_x + arrow_length_ * cos(angle + arrow_degrees_);
	double y2 = end_y + arrow_length_ * sin(angle + arrow_degrees_);
	double w1 = 1.0, w2 = 1.0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_width(cr, fmin(w1, w2) / 2.0);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, start_x, start_y);
	cairo_line_to(cr, (x1 + x2) / 2, (y1 + y2) / 2);
	cairo_stroke(cr);
	cairo_set_line_width(cr, fmin(w1, w2) / 2);
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER);
	cairo_move_to(cr, x1, y1);
	cairo_line_to(cr, x2, y2);
	cairo_line_to(cr, end_x, end_y);
	cairo_line_to(cr, x1, y1);
	cairo_close_path(cr);
	cairo_fill(cr);
}

void draw_ylab(cairo_t *cr, const char *lab, double x, double canvas_height)
{
	cairo_text_extents_t ext;
	// Measure text dimensions
	cairo_text_extents(cr, lab, &ext);
	// Calculate the center of the canvas for the y-axis
	double y_center = canvas_height / 2.0;
	// Adjust position for rotation and centering
	double x_pos = x - (ext.height / 2.0); // Offset to center the rotated text
	double y_pos = y_center + (ext.width / 2.0); // Offset to vertically center the text
	// Apply rotation for vertical text
	cairo_save(cr);
	cairo_translate(cr, x_pos, y_pos);
	cairo_rotate(cr, -M_PI / 2.0); // Rotate 90 degrees counter-clockwise
	// Draw the text
	cairo_move_to(cr, 0, 0);
	cairo_show_text(cr, lab);
	cairo_restore(cr);
}

void draw_yticks(cairo_t *cr, const int ymax)
{
	double x, y;
	double w1 = 1.0, w2 = 1.0;
	int i, j, k, ylm = 0;
	cairo_device_to_user_distance(cr, &w1, &w2);
	cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
	cairo_text_extents_t ext;
	const double dashes[] = {0.75, 5.0, 0.75, 5.0};
	int ndash = sizeof(dashes) / sizeof(dashes[0]);
	char buf[64], buf_ts[64];
	double h = ceil(log10(ymax));
	cairo_text_extents(cr, "m", &ext);
	double x_offset = ext.width;
	for (i = 0; i <= h; ++i)
	{
		memset(buf, 0, 64);
		memset(buf_ts, 0, 64);
		snprintf(buf, 64, "%d", (int)pow(10, h - i));
		int bufl = strlen(buf);
		if (bufl >= 4)
		{
			for (k = 0, j = bufl - 1; j >= 0; --j)
			{
				buf_ts[k++] = buf[j];
				if ((bufl - 1 - j) % 3 == 2 && j > 0)
					buf_ts[k++] = ',';
			}
			reverse(buf_ts);
		}
		else
			strncpy(buf_ts, buf, 64);
		// add thousand separator
		cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
		cairo_text_extents(cr, buf_ts, &ext);
		ylm = fmax(ylm, ext.width);
		x = -ext.width - x_offset / 2.5;
		y = 1 - (double)i / (h + 1);
		cairo_move_to(cr, x, DIM_Y - y * DIM_Y + ext.height / 2);
		cairo_show_text(cr, buf_ts);
		// major ticks
		cairo_set_line_width(cr, fmin(w1, w2));
		cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
		cairo_line_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
		cairo_stroke(cr);
		cairo_set_line_width(cr, fmin(w1, w2) / 4);
		cairo_set_dash(cr, dashes, ndash, 0);
		cairo_move_to(cr, x_offset * .75, DIM_Y - y * DIM_Y);
		cairo_line_to(cr, DIM_X, DIM_Y - y * DIM_Y);
		cairo_stroke(cr);
		cairo_set_dash(cr, dashes, 0, 0);
		cairo_set_line_width(cr, fmin(w1, w2) / 2);
		// minor ticks
		for (j = 2; j <= 9 && i < h; ++j)
		{
			cairo_set_line_width(cr, fmin(w1, w2) / 1.5);
			y = (log10((11 - j) * pow(10, i)) + 1)  / (h + 1);
			cairo_move_to(cr, 0, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, x_offset * .375, DIM_Y - y * DIM_Y);
			cairo_stroke(cr);
			cairo_set_line_width(cr, fmin(w1, w2) / 5);
			cairo_set_dash(cr, dashes, ndash, 0);
			cairo_move_to(cr, x_offset * .375, DIM_Y - y * DIM_Y);
			cairo_line_to(cr, DIM_X, DIM_Y - y * DIM_Y);
			cairo_stroke(cr);
			cairo_set_dash(cr, dashes, 0, 0);
			cairo_set_line_width(cr, fmin(w1, w2) / 2);
		}
	}
	cairo_stroke(cr);
	cairo_select_font_face(cr, "Open Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	draw_ylab(cr, "Depth", fmin(-MARGIN / 3.0, -ylm / 2 * 1.25), DIM_Y);
}

void draw_legend(cairo_t *cr)
{
	cairo_set_font_size(cr, 17.5);
	cairo_text_extents_t ext;
	cairo_text_extents(cr, "m", &ext);
	double yshift = -ext.width / 2;
	cairo_set_line_width(cr, 0.5);
	cairo_set_source_rgb(cr, 36 / 255.0, 123 / 255.0, 160 / 255.0);
	cairo_rectangle(cr, 12.5, -12.5 + yshift, 12.5, 12.5);
	cairo_stroke_preserve(cr);
	cairo_set_source_rgba(cr, 36 / 255.0, 123 / 255.0, 160 / 255.0, 0.9);
	cairo_fill(cr);
	cairo_move_to(cr, 27.5, yshift);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_show_text(cr, "w/o dup");
	cairo_set_source_rgb(cr, 112 / 255.0, 193 / 255.0, 179 / 255.0);
	cairo_rectangle(cr, 107.5, -12.5 + yshift, 12.5, 12.5);
	cairo_stroke_preserve(cr);
	cairo_set_source_rgba(cr, 112 / 255.0, 193 / 255.0, 179 / 255.0, 0.9);
	cairo_fill(cr);
	cairo_move_to(cr, 122.5, yshift);
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_show_text(cr, "w/ dup");
}
