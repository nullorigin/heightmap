/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/functions.hpp"

namespace hmap
{

//----------------------------------------------------------------------
// base class
//----------------------------------------------------------------------

HMAP_FCT_XY_TYPE Function::get_delegate() const
{
  return this->delegate;
}

float Function::get_value(float x, float y, float ctrl_param) const
{
  return this->delegate(x, y, ctrl_param);
}

void Function::set_delegate(HMAP_FCT_XY_TYPE new_delegate)
{
  this->delegate = std::move(new_delegate);
}

//----------------------------------------------------------------------
// derived from Function class
//----------------------------------------------------------------------

ArrayFunction::ArrayFunction(hmap::Array array, Vec2<float> kw, bool periodic)
    : Function(), kw(kw), array(array)
{
  if (periodic)
    this->set_delegate(
        [this](float x, float y, float)
        {
          float xp = 0.5f * this->kw.x * x;
          float yp = 0.5f * this->kw.y * y;

          xp = 2.f * std::abs(xp - int(xp));
          yp = 2.f * std::abs(yp - int(yp));

          xp = xp < 1.f ? xp : 2.f - xp;
          yp = yp < 1.f ? yp : 2.f - yp;

          smoothstep5(xp);
          smoothstep5(yp);

          float xg = xp * (this->array.shape.x - 1);
          float yg = yp * (this->array.shape.y - 1);
          int   i = (int)xg;
          int   j = (int)yg;
          return this->array.get_value_bilinear_at(i, j, xg - i, yg - j);
        });
  else
    this->set_delegate(
        [this](float x, float y, float)
        {
          float xp = std::clamp(this->kw.x * x,
                                0.f,
                                1.f - std::numeric_limits<float>::min());
          float yp = std::clamp(this->kw.y * y,
                                0.f,
                                1.f - std::numeric_limits<float>::min());

          xp = xp - int(xp);
          yp = yp - int(yp);

          float xg = xp * (this->array.shape.x - 1);
          float yg = yp * (this->array.shape.y - 1);
          int   i = (int)xg;
          int   j = (int)yg;
          return this->array.get_value_bilinear_at(i, j, xg - i, yg - j);
        });
}

BumpFunction::BumpFunction(float gain, Vec2<float> center)
    : Function(), center(center)
{
  this->set_gain(gain);
  this->set_delegate(
      [this](float x, float y, float)
      {
        float dx = x - this->center.x;
        float dy = y - this->center.y;

        float r2 = dx * dx + dy * dy;
        return r2 > 0.25f ? 0.f
                          : std::pow(std::exp(-1.f / (1.f - 4.f * r2)),
                                     this->inv_gain);
      });
}

WaveDuneFunction::WaveDuneFunction(Vec2<float> kw,
                                   float       angle,
                                   float       xtop,
                                   float       xbottom,
                                   float       phase_shift)
    : Function(), kw(kw), xtop(xtop), xbottom(xbottom), phase_shift(phase_shift)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * x + sa * this->kw.y * y;
        float xp = std::fmod(r + this->phase_shift + 10.f, 1.f);
        float yp = 0.f;

        if (xp < this->xtop)
        {
          float r = xp / this->xtop;
          yp = r * r * (3.f - 2.f * r);
        }
        else if (xp < this->xbottom)
        {
          float r = (xp - this->xbottom) / (this->xtop - this->xbottom);
          yp = r * r * (2.f - r);
        }
        return yp;
      });
}

WaveSineFunction::WaveSineFunction(Vec2<float> kw,
                                   float       angle,
                                   float       phase_shift)
    : Function(), kw(kw), phase_shift(phase_shift)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * x + sa * this->kw.y * y;
        return std::cos(2.f * M_PI * r + this->phase_shift);
      });
}

WaveSquareFunction::WaveSquareFunction(Vec2<float> kw,
                                       float       angle,
                                       float       phase_shift)
    : Function(), kw(kw), phase_shift(phase_shift)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * x + sa * this->kw.y * y;
        return r = 2.f * (int)r - (int)(2.f * r) + 1.f;
      });
}

WaveTriangularFunction::WaveTriangularFunction(Vec2<float> kw,
                                               float       angle,
                                               float       slant_ratio,
                                               float       phase_shift)
    : Function(), kw(kw), slant_ratio(slant_ratio), phase_shift(phase_shift)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * x + sa * this->kw.y * y;

        r = r - (int)r;
        if (r < this->slant_ratio)
          r /= this->slant_ratio;
        else
          r = 1.f - (r - this->slant_ratio) / (1.f - this->slant_ratio);
        return r * r * (3.f - 2.f * r);
      });
}

} // namespace hmap