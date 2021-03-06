/* lorina: C++ parsing library
 * Copyright (C) 2018  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file blif.hpp
  \brief Implements blif parser

  \author Heinz Riener
*/

#pragma once

#include <lorina/common.hpp>
#include <lorina/diagnostics.hpp>
#include <lorina/detail/utils.hpp>
#include <regex>
#include <iostream>

namespace lorina
{

/*! \brief A reader visitor for the BLIF format.
 *
 * Callbacks for the BLIF (Berkeley Logic Interchange Format) format.
 */
class blif_reader
{
public:
  /*! \brief Type of the output cover as truth table. */
  using output_cover_t = std::vector<std::pair<std::string, std::string>>;

public:
  /*! \brief Callback method for parsed model.
   *
   * \param model_name Name of the model
   */
  virtual void on_model( const std::string& model_name ) const
  {
    (void)model_name;
  }

  /*! \brief Callback method for parsed input.
   *
   * \param name Input name
   */
  virtual void on_input( const std::string& name ) const
  {
    (void)name;
  }

  /*! \brief Callback method for parsed output.
   *
   * \param name Output name
   */
  virtual void on_output( const std::string& name ) const
  {
    (void)name;
  }

  /*! \brief Callback method for parsed gate.
   *
   * \param inputs A list of input names
   * \param output Name of output of the gate
   * \param cover N-input, 1-output PLA description of the logic function corresponding to the logic gate
   */
  virtual void on_gate( const std::vector<std::string>& inputs, const std::string& output, const output_cover_t& cover ) const
  {
    (void)inputs;
    (void)output;
    (void)cover;
  }

  /*! \brief Callback method for parsed end.
   *
   */
  virtual void on_end() const {}

  /*! \brief Callback method for parsed comment.
   *
   * \param comment Comment
   */
  virtual void on_comment( const std::string& comment ) const
  {
    (void)comment;
  }
}; /* blif_reader */

/*! \brief A BLIF reader for prettyprinting BLIF.
 *
 * Callbacks for prettyprinting of BLIF.
 *
 */
class blif_pretty_printer : public blif_reader
{
public:
  /*! \brief Constructor of the BLIF pretty printer.
   *
   * \param os Output stream
   */
  blif_pretty_printer( std::ostream& os = std::cout )
      : _os( os )
  {
  }

  virtual void on_model( const std::string& model_name ) const
  {
    _os << ".model " << model_name << std::endl;
  }

  virtual void on_input( const std::string& name ) const
  {
    if ( first_input ) _os << std::endl << ".inputs ";
    _os << name << " ";
    first_input = false;
  }

  virtual void on_output( const std::string& name ) const
  {
    if ( first_output ) _os << std::endl << ".output ";
    _os << name << " ";
    first_output = false;
  }

  virtual void on_gate( const std::vector<std::string>& inputs, const std::string& output, const output_cover_t& cover ) const
  {
    _os << std::endl << fmt::format( ".names {0} {1}", detail::join( inputs, "," ), output ) << std::endl;
    for ( const auto& c : cover )
    {
      _os << c.first << ' ' << c.second << std::endl;
    }
  }

  virtual void on_end() const
  {
    _os << std::endl << ".end" << std::endl;
  }

  virtual void on_comment( const std::string& comment ) const
  {
    _os << std::endl << "# " << comment << std::endl;
  }

  mutable bool first_input = true; /*!< Predicate that is true until the first input was parsed */
  mutable bool first_output = true; /*!< Predicate that is true until the first output was parsed */
  std::ostream& _os; /*!< Output stream */
}; /* blif_pretty_printer */

namespace blif_regex
{
static std::regex model( R"(.model\s+(.*))" );
static std::regex inputs( R"(.inputs\s+(.*))" );
static std::regex outputs( R"(.outputs\s+(.*))" );
static std::regex names( R"(.names\s+(.*))" );
static std::regex line_of_truthtable( R"(([01\-]*)\s*([01\-]))" );
static std::regex end( R"(.end)" );
} // namespace blif_regex

/*! \brief Reader function for the BLIF format.
 *
 * Reads BLIF format from a stream and invokes a callback
 * method for each parsed primitive and each detected parse error.
 *
 * \param in Input stream
 * \param reader A BLIF reader with callback methods invoked for parsed primitives
 * \param diag An optional diagnostic engine with callback methods for parse errors
 * \return Success if parsing have been successful, or parse error if parsing have failed
 */
inline return_code read_blif( std::istream& in, const blif_reader& reader, diagnostic_engine* diag = nullptr )
{
  return_code result = return_code::success;

  std::smatch m;
  detail::foreach_line_in_file_escape( in, [&]( std::string line ) {
    redo:
      /* empty line or comment */
      if ( line.empty() || line[0] == '#' )
        return true;

      /* names */
      if ( std::regex_search( line, m, blif_regex::names ) )
      {
        auto args = detail::split( detail::trim_copy( m[1] ), " " );
        const auto output = *( args.rbegin() );
        args.pop_back();

        std::vector<std::pair<std::string, std::string>> tt;
        std::string copy_line;
        detail::foreach_line_in_file_escape( in, [&]( const std::string& line ) {
          copy_line = line;

          /* terminate when the next '.' declaration starts */
          if ( line[0] == '.' )
            return false;

          /* empty line or comment */
          if ( line.empty() || line[0] == '#' )
            return true;

          if ( std::regex_search( line, m, blif_regex::line_of_truthtable ) )
          {
            tt.push_back( std::make_pair( detail::trim_copy( m[1] ), detail::trim_copy( m[2] ) ) );
            return true;
          }

          return false;
        } );

        reader.on_gate( args, output, tt );

        if ( in.eof() )
        {
          return true;
        }
        else
        {
          /* restart parsing with the line of the subparser */
          line = copy_line;
          goto redo;
        }
      }

      /* .model <string> */
      if ( std::regex_search( line, m, blif_regex::model ) )
      {
        reader.on_model( detail::trim_copy( m[1] ) );
        return true;
      }

      /* .inputs <list of whitespace separated strings> */
      if ( std::regex_search( line, m, blif_regex::inputs ) )
      {
        for ( const auto& input : detail::split( detail::trim_copy( m[1] ), " " ) )
        {
          reader.on_input( detail::trim_copy( input ) );
        }
        return true;
      }

      /* .outputs <list of whitespace separated strings> */
      if ( std::regex_search( line, m, blif_regex::outputs ) )
      {
        for ( const auto& output : detail::split( detail::trim_copy( m[1] ), " " ) )
        {
          reader.on_output( detail::trim_copy( output ) );
        }
        return true;
      }

      /* .end */
      if ( std::regex_search( line, m, blif_regex::end ) )
      {
        reader.on_end();
        return true;
      }

      if ( diag )
      {
        diag->report( diagnostic_level::error,
                      fmt::format( "cannot parse line `{0}`", line ) );
      }

      result = return_code::parse_error;
      return true;
  } );

  return result;
}

/*! \brief Reader function for BLIF format.
 *
 * Reads binary BLIF format from a file and invokes a callback
 * method for each parsed primitive and each detected parse error.
 *
 * \param filename Name of the file
 * \param reader A BLIF reader with callback methods invoked for parsed primitives
 * \param diag An optional diagnostic engine with callback methods for parse errors
 * \return Success if parsing have been successful, or parse error if parsing have failed
 */
inline return_code read_blif( const std::string& filename, const blif_reader& reader, diagnostic_engine* diag = nullptr )
{
  std::ifstream in( detail::word_exp_filename( filename ), std::ifstream::in );
  return read_blif( in, reader, diag );
}

} // namespace lorina
