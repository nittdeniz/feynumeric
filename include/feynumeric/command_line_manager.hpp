#ifndef FEYNUMERIC_COMMAND_LINE_MANAGER_HPP
#define FEYNUMERIC_COMMAND_LINE_MANAGER_HPP


#include <map>
#include <string>
#include <vector>

namespace Feynumeric
{

	struct Command
	{
		std::string cmd;
		std::string default_value;
		std::string description;
		bool is_mandatory;
	};

	class Command_Line_Manager
	{
	private:
		std::vector<Command> commands;
		std::map<std::string, std::string> arguments;
	public:
		Command_Line_Manager(int argc, char** argv);

		void register_command(const std::string& cmd, bool is_mandatory= false, const std::string& description = "");
		void register_command(std::string const& cmd, std::string const& default_value, std::string const& description = "");

		float as_float(const std::string& arg) const;

		double as_double(std::string const& arg) const;

		int as_int(const std::string& arg) const;

		std::string as_string(const std::string& arg) const;

		bool is_enabled(const std::string& option);

		void crash_on_missing_mandatory_command() const;
	};
}

#endif //FEYNUMERIC_COMMAND_LINE_MANAGER_HPP
