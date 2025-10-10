"""Tools for checking the status of calculations in a calculation directory from a
*status.json* file."""

import os
import pathlib
import typing

from .json_io import read_optional, safe_dump


class CalcDir:
    """Class for checking the status of calculations in a calculation directory from a
    *status.json* file.

    Notes
    -----

    Standard *status.json* attributes are:

    - *"status"* - *str*, Calculation status. Standard values are:
      - *"none"* - Calculation not setup or status unknown.
      - *"setup"* - Calculation input files are set up and ready for calculation.
      - *"submitted"* - Calculation job submitted.
      - *"started"* - Job has begun running.
      - *"canceled"* - Job has been canceled, before or after starting.
      - *"stopped"* - Job has been stopped, after it began running, and is not complete.
        This may be due to an error, hitting the walltime, a canceled job has stopped,
        or some other reason.
      - *"complete"* - Job has completed successfully.
    - *"jobid"* - *str*, The submitted job ID.
    - *"starttime"* - *str*, The time the job started, in ISO 8601 format
      (YYYY-MM-DDTHH:mm:ss)
    - *"stoptime"* - *str*, The time the job was stopped, successfully or not, in
      ISO 8601 format.

    """

    def __init__(
        self,
        path: typing.Optional[pathlib.Path] = None,
        name_parts: int = 2,
    ):
        self.path = path
        """pathlib.Path: The calculation directory path. May be relative or absolute."""

        self.name_parts = name_parts
        """int: The number of parts of the calculation directory path to use for the
        name. 
        
        For example, if the path is `/path/to/calc/dir` and `name_parts` is 2, the
        name will be `calc/dir`. If `name_parts` is 3, the name will be
        `to/calc/dir`. If `name_parts` is greater than the number of parts in the path,
        the entire path is used.
        """

        # Cache for reading calc status.json file
        self._status_data: typing.Optional[dict] = None
        self._status_data_mtime: typing.Optional[float] = None

    @property
    def name(self):
        """str: A name for the calculation directory, based on the last `name_parts`
        parts of the path."""
        parts = self.path.parts
        _name_parts = min(self.name_parts, len(parts))
        return os.path.join(*parts[-_name_parts:])

    @property
    def status_data(self) -> typing.Optional[dict]:
        """Optional[dict]: Contents of the `status.json` file in the
        calculation directory, if it exists.

        If the file contents have changed since the last time it was read, the file
        is read again. If no calculation directory is set or the file does not exist,
        returns None.
        """
        if self.path is None:
            # If no calculation directory is set, return None
            self._status_data = None
            self._status_data_mtime = None
            return None
        status_path = self.path / "status.json"
        if not status_path.exists():
            self._status_data = None
            self._status_data_mtime = None
            return None
        curr_mtime = status_path.stat().st_mtime
        if self._status_data_mtime != curr_mtime:
            # If the mtime has changed, read the file again
            self._status_data = read_optional(status_path)
            self._status_data_mtime = curr_mtime
        return self._status_data

    def commit_status(
        self,
        quiet: bool = False,
    ):
        """Write the current status data to the `status.json` file in the calculation
        directory.

        Parameters
        ----------
        quiet: bool = False
            If True, do not print a message indicating that the file is being written.

        """
        safe_dump(self._status_data, self.path / "status.json", force=True, quiet=quiet)

    @property
    def status(self) -> str:
        """str: The status of the configuration's calculation, as determined by
        checking a `status.json` file in the configuration's calculation directory.

        Returns the value of the `"status"` attribute in the `status.json` file. If
        no `clex` is given for the parent selection or the `status.json` file does
        not exist in the calculation directory, returns "none". If the file does
        exist, but the `"status"` attribute is not present, is not a str, or otherwise
        can't be read, an exception is raised.
        """
        data = self.status_data
        if data is None:
            return "none"
        status = data.get("status")
        if status is None:
            return "none"
        elif not isinstance(status, str):
            raise ValueError(f"Invalid status: expected str, got {type(status)}")
        return status

    @property
    def jobid(self) -> str:
        """str: The calculation's job ID, as determined by checking a `status.json`
        file in the configuration's calculation directory.

        Returns the value of the `"jobid"` attribute in the `status.json` file. If
        no `clex` is given for the parent selection or the `status.json` file does
        not exist in the calculation directory, or a `"jobid"` is not present,
        returns "none". If the file does exist, but the `"status"` attribute is not
        present, is not a str, or otherwise can't be read, an exception is raised.
        """
        data = self.status_data
        if data is None:
            return "none"
        jobid = data.get("jobid")
        if jobid is None:
            return "none"
        elif not isinstance(jobid, str):
            raise ValueError(f"Invalid jobid: expected str, got {type(jobid)}")
        return jobid

    @property
    def starttime(self) -> str:
        """Optional[str]: The calculation's start time, as determined by checking a
        `status.json` file in the configuration's calculation directory.

        Returns the value of the `"starttime"` attribute in the `status.json` file.
        If no `clex` is given for the parent selection or the `status.json` file does
        not exist in the calculation directory, returns "none". If the file does exist,
        but the `"starttime"` attribute is not present, is not a str, or otherwise
        can't be read, an exception is raised.
        """
        data = self.status_data
        if data is None:
            return "none"
        starttime = data.get("starttime")
        if starttime is None:
            return "none"
        elif not isinstance(starttime, str):
            raise ValueError(f"Invalid starttime: expected str, got {type(starttime)}")
        return starttime

    @property
    def stoptime(self) -> str:
        """Optional[str]: The calculation's stop time, as determined by checking a
        `status.json` file in the configuration's calculation directory.

        Returns the value of the `"stoptime"` attribute in the `status.json` file.
        If no `clex` is given for the parent selection or the `status.json` file does
        not exist in the calculation directory, returns "none". If the file does exist,
        but the `"stoptime"` attribute is not present, is not a str, or otherwise
        can't be read, an exception is raised.
        """
        data = self.status_data
        if data is None:
            return "none"
        stoptime = data.get("stoptime")
        if stoptime is None:
            return "none"
        elif not isinstance(stoptime, str):
            raise ValueError(f"Invalid stoptime: expected str, got {type(stoptime)}")
        return stoptime

    @property
    def runtime(self) -> str:
        """str: The calculation's runtime in HH:MM:SS format, as determined by
        checking a `status.json` file in the configuration's calculation directory.

        Calculates the runtime as the difference between `stoptime` and
        `starttime`, if both are available, and formats the value as a string. If
        either is not available, returns "none".
        """
        starttime = self.starttime
        stoptime = self.stoptime
        if starttime == "none":
            return "none"

        # If stoptime is "none", we assume the calculation is still running and
        # set stoptime to the current time.
        running = False
        if stoptime == "none":
            from datetime import datetime

            stoptime = datetime.now().isoformat()
            running = True

        # starttime and stoptime are expected to be in format generated by
        # `$(date +%Y-%m-%dT%H:%M:%S)` in a bash script, which is ISO 8601 format.
        from datetime import datetime, timedelta

        try:
            start_dt = datetime.fromisoformat(starttime)
            stop_dt = datetime.fromisoformat(stoptime)
        except ValueError as e:
            raise ValueError(
                f"Invalid datetime format in status.json: "
                f"starttime='{starttime}', stoptime='{stoptime}'"
            ) from e
        runtime = stop_dt - start_dt

        # Runtime is in timedelta format. Convert to D-HH:MM:SS format:
        def formatted(runtime: timedelta) -> str:
            total_seconds = int(runtime.total_seconds())

            days = total_seconds // (24 * 3600)
            total_seconds %= 24 * 3600
            hours = total_seconds // 3600
            total_seconds %= 3600
            minutes = total_seconds // 60
            seconds = total_seconds % 60

            if days > 0:
                return f"{days}-{hours:02}:{minutes:02}:{seconds:02}"
            else:
                if hours > 0:
                    return f"{hours}:{minutes:02}:{seconds:02}"
                elif minutes > 0:
                    return f"{minutes}:{seconds:02}"
                else:
                    # Handles cases with only seconds, e.g., 30 seconds -> 0:30
                    return f"0:{seconds:02}"

        formatted_time = formatted(runtime)
        if running:
            formatted_time = f"{formatted_time}+"
        return formatted_time
