from setuptools import setup

entry_points = {
    "console_scripts": [
        "fss_check_daily_fss=fss_check.scripts.daily_fss:main",
    ]
}

setup(
    name="fss_check",
    author="Tom Aldcroft",
    description="Fine Sun Sensor monitoring and trending",
    author_email="taldcroft@cfa.harvard.edu",
    use_scm_version=True,
    setup_requires=["setuptools_scm", "setuptools_scm_git_archive"],
    zip_safe=False,
    entry_points=entry_points,
    packages=["fss_check", "fss_check.scripts"],
    package_data={
        "fss_check": ["fss_check_config.yml", "data/index.html", "task_schedule.cfg"]
    },
)
