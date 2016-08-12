/*
* Copyright (C) 2016 Swift Navigation Inc.
* Contact: Pasi Miettinen <pasi.miettinen@exafore.com>
*
* This source is subject to the license found in the file 'LICENSE'
* which must be be distributed together with this source. All other
* rights reserved.
*
* THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
* KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
* PURPOSE.
*/

#ifndef NOVATELPYTHON_H
#define NOVATELPYTHON_H

#include "novatel.h"

#include <boost/python.hpp>
using namespace boost::python;

namespace novatel {

#define SIZEOF_ARRAY(a) (sizeof(a) / sizeof(a[0]))

class NovatelPython : public Novatel
{
public:
  NovatelPython(std::string id);
  NovatelPython(std::string id, std::string raw_log_name);
  ~NovatelPython();

  template<typename T>
  static boost::python::tuple arr2tuple(T *arr, uint16_t len) {
    boost::python::list ret;
    for (uint8_t i = 0; i < len; i++) {
      ret.append(arr[i]);
    }
    return boost::python::tuple(ret);
  }

  PyObject *
    python_set_best_position_callback(PyObject *callback);
  PyObject *
    python_set_best_position_ecef_callback(PyObject *callback);
  PyObject *
    python_set_range_measurements_callback(PyObject *callback);
  PyObject *
    python_set_raw_gps_word_callback(PyObject *callback);
  PyObject *
    python_set_gps_ephem_callback(PyObject *callback);
  PyObject *
    python_set_gps_almanac_callback(PyObject *callback);
  PyObject *
    python_set_best_sats_callback(PyObject *callback);

  void best_position_callback(Position const &best_pos,
    double const &timestamp);
  void best_position_ecef_callback(PositionEcef const &best_pos_ecef,
    double const &timestamp);
  void range_measurements_callback(RangeMeasurements const &range_meas,
    double const &timestamp);
  void raw_gps_word_callback(RawGpsWord const &raw_gps_word,
    double const &timestamp);
  void gps_ephem_callback(GpsEphemeris const &gps_ephem,
    double const &timestamp);
  void gps_almanac_callback(Almanac const &gps_almanac,
    double const &timestamp);
  void best_sats_callback(BestSats const &best_sats,
    double const &timestamp);

private:
  std::string m_id;

  PyObject *python_best_position_callback_;
  PyObject *python_best_position_ecef_callback_;
  PyObject *python_range_measurements_callback_;
  PyObject *python_raw_gps_word_callback_;
  PyObject *python_gps_ephem_callback_;
  PyObject *python_gps_almanac_callback_;
  PyObject *python_best_sats_callback_;

  template<typename T>
  void call_python(PyObject *callable, T const &data,
    double const &timestamp);
};

class Oem4BinaryHeaderPython {
public:
  Oem4BinaryHeaderPython(Oem4BinaryHeader const &hdr) {m_hdr = hdr;}
  uint8_t sync1() {return m_hdr.sync1;}
  uint8_t sync2() {return m_hdr.sync2;}
  uint8_t sync3() {return m_hdr.sync3;}
  uint8_t header_length() {return m_hdr.header_length;}
  uint16_t message_id() {return m_hdr.message_id;}
  uint8_t message_type() {
    uint8_t ret = m_hdr.message_type.reserved << 3;
    ret |= m_hdr.message_type.format << 1;
    ret |= m_hdr.message_type.response;
    return ret;
  }
  uint8_t port_address() {return m_hdr.port_address;}
  uint16_t message_length() {return m_hdr.message_length;}
  uint16_t sequence() {return m_hdr.sequence;}
  uint8_t idle() {return m_hdr.idle;}
  uint8_t time_status() {return m_hdr.time_status;}
  uint16_t gps_week() {return m_hdr.gps_week;}
  uint32_t gps_millisecs() {return m_hdr.gps_millisecs;}
  uint32_t status() {return m_hdr.status;}
  uint16_t Reserved() {return m_hdr.Reserved;}
  uint16_t version() {return m_hdr.version;}
private:
  Oem4BinaryHeader m_hdr;
};

class PositionPython {
public:
  PositionPython(Position const &position) {
    m_pos = position;
  }
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_pos.header);
  }
  int32_t solution_status() {return m_pos.solution_status;}
  int32_t position_type() {return m_pos.position_type;}
  double latitude() {return m_pos.latitude;}
  double longitude() {return m_pos.longitude;}
  double height() {return m_pos.height;}
  float undulation() {return m_pos.undulation;}
  int32_t datum_id() {return m_pos.datum_id;}
  float latitude_standard_deviation() {
    return m_pos.latitude_standard_deviation;
  }
  float longitude_standard_deviation() {
    return m_pos.longitude_standard_deviation;
  }
  float height_standard_deviation() {
    return m_pos.height_standard_deviation;
  }
  boost::python::tuple base_station_id() {
    return NovatelPython::arr2tuple(m_pos.base_station_id,
      sizeof(m_pos.base_station_id));
  }
  float differential_age() {return m_pos.differential_age;}
  float solution_age() {return m_pos.solution_age;}
  uint8_t number_of_satellites() {return m_pos.number_of_satellites;}
  uint8_t number_of_satellites_in_solution() {
    return m_pos.number_of_satellites_in_solution;
  }
  uint8_t num_gps_plus_glonass_l1() {
    return m_pos.num_gps_plus_glonass_l1;
  }
  uint8_t num_gps_plus_glonass_l2() {
    return m_pos.num_gps_plus_glonass_l2;
  }
  uint8_t reserved() {return m_pos.reserved;}
  uint8_t extended_solution_status() {
    return m_pos.extended_solution_status;
  }
  uint8_t reserved2() {return m_pos.reserved2;}
  uint8_t signals_used_mask() {return m_pos.signals_used_mask;}
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_pos.crc, sizeof(m_pos.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    d["time_status"] = m_pos.header.time_status;
    d["gps_week"] = m_pos.header.gps_week;
    d["gps_millisecs"] = m_pos.header.gps_millisecs;
    d["status"] = m_pos.header.status;
    d["solution_status"] = solution_status();
    d["position_type"] = position_type();
    d["latitude"] = latitude();
    d["longitude"] = longitude();
    d["height"] = height();
    d["undulation"] = undulation();
    d["datum_id"] = datum_id();
    d["latitude_standard_deviation"] = latitude_standard_deviation();
    d["longitude_standard_deviation"] = longitude_standard_deviation();
    d["height_standard_deviation"] = height_standard_deviation();
    std::string bsid = (char*)m_pos.base_station_id;
    d["base_station_id"] = bsid;
    d["differential_age"] = differential_age();
    d["solution_age"] = solution_age();
    d["number_of_satellites"] = number_of_satellites();
    d["number_of_satellites_in_solution"] =
      number_of_satellites_in_solution();
    d["num_gps_plus_glonass_l1"] = num_gps_plus_glonass_l1();
    d["num_gps_plus_glonass_l2"] = num_gps_plus_glonass_l2();
    d["reserved"] = reserved();
    d["extended_solution_status"] = extended_solution_status();
    d["reserved2"] = reserved2();
    d["signals_used_mask"] = signals_used_mask();
    return d;
  }
private:
  Position m_pos;
};

class PositionEcefPython {
public:
  PositionEcefPython(PositionEcef const &position) {
    m_pos = position;
  }
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_pos.header);
  }
  int32_t solution_status() {return m_pos.solution_status;}
  int32_t position_type() {return m_pos.position_type;}
  double x_position() {return m_pos.x_position;}
  double y_position() {return m_pos.y_position;}
  double z_position() {return m_pos.z_position;}
  float x_standard_deviation() {return m_pos.x_standard_deviation;}
  float y_standard_deviation() {return m_pos.y_standard_deviation;}
  float z_standard_deviation() {return m_pos.z_standard_deviation;}
  int32_t velocity_status() {return m_pos.velocity_status;}
  int32_t velocity_type() {return m_pos.velocity_type;}
  double x_velocity() {return m_pos.x_velocity;}
  double y_velocity() {return m_pos.y_velocity;}
  double z_velocity() {return m_pos.z_velocity;}
  float x_velocity_standard_deviation() {
    return m_pos.x_velocity_standard_deviation;
  }
  float y_velocity_standard_deviation() {
    return m_pos.y_velocity_standard_deviation;
  }
  float z_velocity_standard_deviation() {
    return m_pos.z_velocity_standard_deviation;
  }
  boost::python::tuple base_station_id() {
    return NovatelPython::arr2tuple(m_pos.base_station_id,
      sizeof(m_pos.base_station_id));
  }
  float velocity_latency() {return m_pos.velocity_latency;}
  float differential_age() {return m_pos.solution_age;}
  float solution_age() {return m_pos.position_type;}
  uint8_t number_of_satellites() {return m_pos.number_of_satellites;}
  uint8_t number_of_satellites_in_solution() {
    return m_pos.number_of_satellites_in_solution;
  }
  boost::python::tuple reserved() {
    return NovatelPython::arr2tuple(m_pos.reserved,
      sizeof(m_pos.reserved));
  }
  uint8_t extended_solution_status() {
    return m_pos.extended_solution_status;
  }
  uint8_t reserved2() {return m_pos.reserved2;}
  uint8_t signals_used_mask() {return m_pos.signals_used_mask;}
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_pos.crc, sizeof(m_pos.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    d["time_status"] = m_pos.header.time_status;
    d["gps_week"] = m_pos.header.gps_week;
    d["gps_millisecs"] = m_pos.header.gps_millisecs;
    d["status"] = m_pos.header.status;
    d["solution_status"] = solution_status();
    d["position_type"] = position_type();
    d["x_position"] = x_position();
    d["y_position"] = y_position();
    d["z_position"] = z_position();
    d["x_standard_deviation"] = x_standard_deviation();
    d["y_standard_deviation"] = y_standard_deviation();
    d["z_standard_deviation"] = z_standard_deviation();
    d["velocity_status"] = velocity_status();
    d["velocity_type"] = velocity_type();
    d["x_velocity_standard_deviation"] =
      x_velocity_standard_deviation();
    d["y_velocity_standard_deviation"] =
      y_velocity_standard_deviation();
    d["z_velocity_standard_deviation"] =
      z_velocity_standard_deviation();
    std::string bsid = (char*)m_pos.base_station_id;
    d["base_station_id"] = bsid;
    d["velocity_latency"] = velocity_latency();
    d["differential_age"] = differential_age();
    d["solution_age"] = solution_age();
    d["number_of_satellites"] = number_of_satellites();
    d["number_of_satellites_in_solution"] =
      number_of_satellites_in_solution();
    d["extended_solution_status"] = extended_solution_status();
    d["signals_used_mask"] = signals_used_mask();
    return d;
  }
private:
  PositionEcef m_pos;
};

class ChannelStatusPython {
public:
  ChannelStatusPython(ChannelStatus const &channel_status) {
    m_status = channel_status;
  }
  uint8_t tracking_state() {return m_status.tracking_state;}
  uint8_t sv_chan_num() {return m_status.sv_chan_num;}
  uint8_t phase_lock_flag() {return m_status.phase_lock_flag;}
  uint8_t parity_known_flag() {return m_status.parity_known_flag;}
  uint8_t code_locked_flag() {return m_status.code_locked_flag;}
  uint8_t correlator_type() {return m_status.correlator_type;}
  uint8_t satellite_sys() {return m_status.satellite_sys;}
  uint8_t reserved1() {return m_status.reserved1;}
  uint8_t grouping() {return m_status.grouping;}
  uint8_t signal_type() {return m_status.signal_type;}
  uint8_t forward_err_correction() {
    return m_status.forward_err_correction;
  }
  uint8_t primary_L1_chan() {return m_status.primary_L1_chan;}
  uint8_t carrier_phase_meas() {return m_status.carrier_phase_meas;}
  uint8_t reserved2() {return m_status.reserved2;}
  uint8_t prn_lock_flag() {return m_status.prn_lock_flag;}
  uint8_t channel_assignment() {return m_status.channel_assignment;}

  uint32_t packed() {
    uint32_t packed = m_status.tracking_state << 27;
    packed |= m_status.sv_chan_num << 22;
    packed |= m_status.phase_lock_flag << 21;
    packed |= m_status.parity_known_flag << 20;
    packed |= m_status.code_locked_flag << 19;
    packed |= m_status.correlator_type << 16;
    packed |= m_status.satellite_sys << 13;
    packed |= m_status.reserved1 << 12;
    packed |= m_status.grouping << 11;
    packed |= m_status.signal_type << 6;
    packed |= m_status.forward_err_correction << 5;
    packed |= m_status.primary_L1_chan << 4;
    packed |= m_status.carrier_phase_meas << 3;
    packed |= m_status.reserved2 << 2;
    packed |= m_status.prn_lock_flag << 1;
    packed |= m_status.channel_assignment;
    return packed;
  }
private:
  ChannelStatus m_status;
};

class RangeDataPython {
public:
  RangeDataPython(RangeData const &data) {m_data = data;}
  
  uint16_t satellite_prn() {return m_data.satellite_prn;}
  uint16_t glonass_frequency() {return m_data.glonass_frequency;}
  double pseudorange() {return m_data.pseudorange;}
  float pseudorange_standard_deviation() {
    return m_data.pseudorange_standard_deviation;
  }
  double accumulated_doppler() {return m_data.accumulated_doppler;}
  float accumulated_doppler_std_deviation() {
    return m_data.accumulated_doppler_std_deviation;
  }
  float doppler() {return m_data.doppler;}
  float carrier_to_noise() {return m_data.carrier_to_noise;}
  float locktime() {return m_data.locktime;}
  uint32_t channel_status_packed() {
    return ChannelStatusPython(m_data.channel_status).packed();
  }
  ChannelStatusPython channel_status() {
    return ChannelStatusPython(m_data.channel_status);
  }
private:
  RangeData m_data;
};

class RangeMeasurementsPython {
public:
  RangeMeasurementsPython(RangeMeasurements const &meas) {
    m_measurements = meas;
  }
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_measurements.header);
  }
  int32_t number_of_observations() {
    return m_measurements.number_of_observations;
  }
  boost::python::tuple range_data() {
    boost::python::list ret;
    for (uint16_t i = 0; i < number_of_observations(); i++) {
      ret.append(RangeDataPython(m_measurements.range_data[i]));
    }
    return boost::python::tuple(ret);
  }
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_measurements.crc,
      sizeof(m_measurements.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    boost::python::list l;
    d["number_of_observations"] = number_of_observations();
    for (uint16_t i = 0; i < number_of_observations(); i++) {
      boost::python::dict meas;
      meas["time_status"] = m_measurements.header.time_status;
      meas["gps_week"] = m_measurements.header.gps_week;
      meas["gps_millisecs"] = m_measurements.header.gps_millisecs;
      meas["status"] = m_measurements.header.status;
      meas["satellite_prn"] =
        m_measurements.range_data[i].satellite_prn;
      meas["glonass_frequency"] =
        m_measurements.range_data[i].glonass_frequency;
      meas["pseudorange"] = m_measurements.range_data[i].pseudorange;
      meas["pseudorange_standard_deviation"] =
        m_measurements.range_data[i].pseudorange_standard_deviation;
      meas["accumulated_doppler"] =
        m_measurements.range_data[i].accumulated_doppler;
      meas["accumulated_doppler_std_deviation"] =
        m_measurements.range_data[i].accumulated_doppler_std_deviation;
      meas["doppler"] = m_measurements.range_data[i].doppler;
      meas["carrier_to_noise"] =
        m_measurements.range_data[i].carrier_to_noise;
      meas["locktime"] = m_measurements.range_data[i].locktime;
      meas["satellite_sys"] =
        m_measurements.range_data[i].channel_status.satellite_sys;
      meas["channel_status_packed"] =
        RangeDataPython(m_measurements.range_data[i]).\
        channel_status_packed();
      l.append(meas);
    }
    d["range_data"] = l;
    return d;
  }
private:
  RangeMeasurements m_measurements;
};

class RawGpsWordPython {
public:
  RawGpsWordPython(RawGpsWord const &word) {m_word = word;}
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_word.header);
  }
  uint32_t prn() {return m_word.prn;}
  boost::python::tuple nav_word() {
    return NovatelPython::arr2tuple(m_word.nav_word,
      sizeof(m_word.nav_word));
  }
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_word.crc, sizeof(m_word.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    boost::python::list l;
    d["time_status"] = m_word.header.time_status;
    d["gps_week"] = m_word.header.gps_week;
    d["gps_millisecs"] = m_word.header.gps_millisecs;
    d["status"] = m_word.header.status;
    d["prn"] = prn();
    std::stringstream ss;
    ss << std::hex << std::uppercase << std::setfill('0');
    for (uint8_t i = 0; i < sizeof(m_word.nav_word); i++) {
      ss << std::setw(2) << (uint8_t)m_word.nav_word[i];
    }
    d["nav_word"] = ss.str();
    return d;
  }
private:
  RawGpsWord m_word;
};

class GpsEphemerisPython {
public:
  GpsEphemerisPython(GpsEphemeris const &ephe) {m_ephe = ephe;}
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_ephe.header);
  }
  uint32_t prn() {return m_ephe.prn;}
  double time_of_week() {return m_ephe.time_of_week;}
  uint32_t health() {return m_ephe.health;}
  uint32_t issue_of_ephemeris_1() {return m_ephe.issue_of_ephemeris_1;}
  uint32_t issue_of_ephemeris_2() {return m_ephe.issue_of_ephemeris_2;}
  uint32_t gps_week() {return m_ephe.gps_week;}
  uint32_t z_count_week() {return m_ephe.z_count_week;}
  double time_of_ephemeris() {return m_ephe.time_of_ephemeris;}
  double semi_major_axis() {return m_ephe.semi_major_axis;}
  double mean_motion_difference() {
    return m_ephe.mean_motion_difference;
  }
  double anomoly_reference_time() {
    return m_ephe.anomoly_reference_time;
  }
  double eccentricity() {return m_ephe.eccentricity;}
  double omega() {return m_ephe.omega;}
  double latitude_cosine() {return m_ephe.latitude_cosine;}
  double latitude_sine() {return m_ephe.latitude_sine;}
  double orbit_radius_cosine() {return m_ephe.orbit_radius_cosine;}
  double orbit_radius_sine() {return m_ephe.orbit_radius_sine;}
  double inclination_cosine() {return m_ephe.inclination_cosine;}
  double inclination_sine() {return m_ephe.inclination_sine;}
  double inclination_angle() {return m_ephe.inclination_angle;}
  double inclination_angle_rate() {
    return m_ephe.inclination_angle_rate;
  }
  double right_ascension() {return m_ephe.right_ascension;}
  double right_ascension_rate() {return m_ephe.right_ascension_rate;}
  uint32_t issue_of_data_clock() {return m_ephe.issue_of_data_clock;}
  double sv_clock_correction() {return m_ephe.sv_clock_correction;}
  double group_delay_difference() {
    return m_ephe.group_delay_difference;
  }
  double clock_aligning_param_0() {
    return m_ephe.clock_aligning_param_0;
  }
  double clock_aligning_param_1() {
    return m_ephe.clock_aligning_param_1;
  }
  double clock_aligning_param_2() {
    return m_ephe.clock_aligning_param_2;
  }
  uint8_t anti_spoofing() {return m_ephe.anti_spoofing;}
  double corrected_mean_motion() {return m_ephe.corrected_mean_motion;}
  double range_accuracy_variance() {
    return m_ephe.range_accuracy_variance;
  }
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_ephe.crc, sizeof(m_ephe.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    d["time_status"] = m_ephe.header.time_status;
    d["gps_week"] = m_ephe.header.gps_week;
    d["gps_millisecs"] = m_ephe.header.gps_millisecs;
    d["status"] = m_ephe.header.status;
    d["prn"] = prn();
    d["time_of_week"] = time_of_week();
    d["health"] = health();
    d["issue_of_ephemeris_1"] = issue_of_ephemeris_1();
    d["issue_of_ephemeris_2"] = issue_of_ephemeris_2();
    d["gps_week"] = gps_week();
    d["z_count_week"] = z_count_week();
    d["time_of_ephemeris"] = time_of_ephemeris();
    d["semi_major_axis"] = semi_major_axis();
    d["mean_motion_difference"] = mean_motion_difference();
    d["anomoly_reference_time"] = anomoly_reference_time();
    d["eccentricity"] = eccentricity();
    d["omega"] = omega();
    d["latitude_cosine"] = latitude_cosine();
    d["latitude_sine"] = latitude_sine();
    d["orbit_radius_cosine"] = orbit_radius_cosine();
    d["orbit_radius_sine"] = orbit_radius_sine();
    d["inclination_cosine"] = inclination_cosine();
    d["inclination_sine"] = inclination_sine();
    d["inclination_angle"] = inclination_angle();
    d["inclination_angle_rate"] = inclination_angle_rate();
    d["right_ascension"] = right_ascension();
    d["right_ascension_rate"] = right_ascension_rate();
    d["issue_of_data_clock"] = issue_of_data_clock();
    d["sv_clock_correction"] = sv_clock_correction();
    d["group_delay_difference"] = group_delay_difference();
    d["clock_aligning_param_0"] = clock_aligning_param_0();
    d["clock_aligning_param_1"] = clock_aligning_param_1();
    d["clock_aligning_param_2"] = clock_aligning_param_2();
    d["anti_spoofing"] = anti_spoofing();
    d["corrected_mean_motion"] = corrected_mean_motion();
    d["range_accuracy_variance"] = range_accuracy_variance();
    return d;
  }
private:
  GpsEphemeris m_ephe;
};

class AlmanacDataPython {
public:
  AlmanacDataPython(AlmanacData const &data) {m_data = data;}
  uint32_t prn() {return m_data.prn;}
	uint32_t ref_week() {return m_data.ref_week;}
	double ref_time() {return m_data.ref_time;}
	double eccentricity() {return m_data.eccentricity;}
	double right_ascension_rate() {return m_data.right_ascension_rate;}
	double right_ascension() {return m_data.right_ascension;}
	double perigee() {return m_data.perigee;}
	double mean_anomoly_of_ref_time() {
    return m_data.mean_anomoly_of_ref_time;
  }
	double clock_aging_param_0() {return m_data.clock_aging_param_0;}
	double clock_aging_param_1() {return m_data.clock_aging_param_1;}
	double corrected_mean_motion() {return m_data.corrected_mean_motion;}
	double semi_major_axis() {return m_data.semi_major_axis;}
	double inclination_angle() {return m_data.inclination_angle;}
	uint32_t sv_configuration() {return m_data.sv_configuration;}
	uint32_t sv_health() {return m_data.sv_health;}
	uint32_t sv_health_from_almanac() {
    return m_data.sv_health_from_almanac;
  }
	uint8_t anti_spoofing() {return m_data.anti_spoofing;}
private:
  AlmanacData m_data;
};

class AlmanacPython {
public:
  AlmanacPython(Almanac const &almanac) {m_almanac = almanac;}
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_almanac.header);
  }
  int32_t number_of_prns() {return m_almanac.number_of_prns;}
  boost::python::tuple data() {
    boost::python::list ret;
    for (uint16_t i = 0; i < number_of_prns(); i++) {
      ret.append(AlmanacDataPython(m_almanac.data[i]));
    }
    return boost::python::tuple(ret);
  }
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_almanac.crc,
      sizeof(m_almanac.crc));
  }
  boost::python::dict as_dict() {
    boost::python::dict d;
    boost::python::list l;
    d["number_of_prns"] = number_of_prns();
    for (uint16_t i = 0; i < number_of_prns(); i++) {
      boost::python::dict alm;
      alm["time_status"] = m_almanac.header.time_status;
      alm["gps_week"] = m_almanac.header.gps_week;
      alm["gps_millisecs"] = m_almanac.header.gps_millisecs;
      alm["status"] = m_almanac.header.status;
      alm["prn"] = m_almanac.data[i].prn;
      alm["ref_week"] = m_almanac.data[i].ref_week;
      alm["ref_time"] = m_almanac.data[i].ref_time;
      alm["eccentricity"] = m_almanac.data[i].eccentricity;
      alm["right_ascension_rate"] = m_almanac.data[i].right_ascension_rate;
      alm["right_ascension"] = m_almanac.data[i].right_ascension;
      alm["perigee"] = m_almanac.data[i].perigee;
      alm["mean_anomoly_of_ref_time"] =
        m_almanac.data[i].mean_anomoly_of_ref_time;
      alm["clock_aging_param_0"] = m_almanac.data[i].clock_aging_param_0;
      alm["clock_aging_param_1"] = m_almanac.data[i].clock_aging_param_1;
      alm["corrected_mean_motion"] = m_almanac.data[i].corrected_mean_motion;
      alm["semi_major_axis"] = m_almanac.data[i].semi_major_axis;
      alm["inclination_angle"] = m_almanac.data[i].inclination_angle;
      alm["sv_configuration"] = m_almanac.data[i].sv_configuration;
      alm["sv_health"] = m_almanac.data[i].sv_health;
      alm["sv_health_from_almanac"] = m_almanac.data[i].sv_health_from_almanac;
      alm["anti_spoofing"] = m_almanac.data[i].anti_spoofing;
      l.append(alm);
    }
    d["almas"] = l;
    return d;
  }
private:
  Almanac m_almanac;
};

class BestSatsPython {
public:
  BestSatsPython(BestSats const &best_sats) {m_sats = best_sats;}
  Oem4BinaryHeaderPython header() {
    return Oem4BinaryHeaderPython(m_sats.header);
  }
  int32_t entries() {return m_sats.entries;}
  boost::python::tuple sats() {
    return NovatelPython::arr2tuple(m_sats.sats, m_sats.entries);
  }
  boost::python::tuple crc() {
    return NovatelPython::arr2tuple(m_sats.crc, sizeof(m_sats.crc));
  }
private:
  BestSats m_sats;
};

}
#endif
