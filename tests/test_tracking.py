import unittest
from unittest.mock import Mock, patch
import logging

from utils.command_data import CommandData
from utils.tracking import execute_command


class TestExecuteCommand(unittest.TestCase):
    def setUp(self):
        self.mock_obj = Mock()
        self.mock_obj.attr1 = 10
        self.mock_obj.attr2 = "initial"
        self.mock_obj.method = Mock(return_value="result")
        self.mock_logger = Mock(spec=logging.Logger)

    @patch("utils.tracking.default_logger")
    def test_simple_tracking(self, mock_default_logger):
        command_data = CommandData(
            obj=self.mock_obj, command="method", tracked_attributes=["attr1"]
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["result"], "result")
        self.assertEqual(result["initial_values"]["attr1"], 10)
        self.assertEqual(result["final_values"]["attr1"], 10)
        self.assertEqual(result["changes"]["attr1"], 0)

    @patch("utils.tracking.default_logger")
    def test_numeric_and_non_numeric_attributes(self, mock_default_logger):
        self.mock_obj.method = Mock(
            side_effect=lambda: setattr(self.mock_obj, "attr1", 15)
        )
        command_data = CommandData(
            obj=self.mock_obj, command="method", tracked_attributes=["attr1", "attr2"]
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["changes"]["attr1"], 5)
        self.assertNotIn("attr2", result["changes"])

    def test_no_tracked_attributes(self):
        command_data = CommandData(
            obj=self.mock_obj, command="method", tracked_attributes=[]
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["result"], "result")
        self.assertEqual(result["initial_values"], {})
        self.assertEqual(result["final_values"], {})
        self.assertEqual(result["changes"], {})

    def test_invalid_attribute(self):
        # Remove the invalid attribute if it exists (it might be auto-created by Mock)
        if hasattr(self.mock_obj, "invalid_attr"):
            delattr(self.mock_obj, "invalid_attr")

        command_data = CommandData(
            obj=self.mock_obj, command="method", tracked_attributes=["invalid_attr"]
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertNotIn("invalid_attr", result["initial_values"])
        self.assertNotIn("invalid_attr", result["final_values"])
        self.assertNotIn("invalid_attr", result["changes"])

    def test_validation_pass(self):
        def validation(obj, initial, final, changes):
            return True

        command_data = CommandData(
            obj=self.mock_obj,
            command="method",
            tracked_attributes=["attr1"],
            validations=[validation],
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["result"], "result")
        self.assertTrue(result["validation_results"]["validation_1"])

    def test_validation_fail(self):
        def validation(obj, initial, final, changes):
            return False

        command_data = CommandData(
            obj=self.mock_obj,
            command="method",
            tracked_attributes=["attr1"],
            validations=[validation],
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["result"], "result")
        self.assertFalse(result["validation_results"]["validation_1"])

    def test_multiple_validations(self):
        def validation1(obj, initial, final, changes):
            return True

        def validation2(obj, initial, final, changes):
            return False

        command_data = CommandData(
            obj=self.mock_obj,
            command="method",
            tracked_attributes=["attr1"],
            validations=[validation1, validation2],
        )
        result = execute_command(command_data, logger=self.mock_logger)
        self.assertEqual(result["result"], "result")
        self.assertTrue(result["validation_results"]["validation_1"])
        self.assertFalse(result["validation_results"]["validation_2"])

    @patch("utils.tracking.default_logger")
    def test_error_handling(self, mock_default_logger):
        self.mock_obj.method = Mock(side_effect=ValueError("Test error"))
        command_data = CommandData(
            obj=self.mock_obj, command="method", tracked_attributes=["attr1"]
        )
        with self.assertRaises(ValueError):
            execute_command(command_data, logger=self.mock_logger)

    def test_callable_command(self):
        def callable_command(obj, *args, **kwargs):
            obj.attr1 += 5
            return "callable result"

        command_data = CommandData(
            obj=self.mock_obj,
            command=callable_command,
            tracked_attributes=["attr1"],
            args=(1, 2),
            kwargs={"key": "value"}
        )
        result = execute_command(command_data, logger=self.mock_logger)
        
        self.assertEqual(result["result"], "callable result")
        self.assertEqual(result["initial_values"]["attr1"], 10)
        self.assertEqual(result["final_values"]["attr1"], 15)
        self.assertEqual(result["changes"]["attr1"], 5)
        self.mock_logger.debug.assert_called_with("Executed command '<function TestExecuteCommand.test_callable_command.<locals>.callable_command at ...>' with result: callable result")


if __name__ == "__main__":
    unittest.main()
