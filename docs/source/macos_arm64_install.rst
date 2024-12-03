Fortran installation in Mac OS (Apple Silicon M1/M2/M3)
-----------------------------------

If you are reading this section, it is because you have an Apple Silicon chip and encountered the following error when trying to access any variable, class or function of the package:

.. code-block:: bash

    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    ImportError: dlopen(/Users/acuna/anaconda3/envs/GASTLI_norosetta/lib/python3.10/site-packages/gastli/dimensions.cpython-310-darwin.so, 0x0002): symbol not found in flat namespace (_f2pyinitdimensions_)

There are two options to fix this error: 1) using an Intel emulation with Rosetta, or 2) using a disk image. We explain both options below.

1) Installation with Rosetta (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of GASTLI on Apple machines with the M1/M2/M3 chip requires Intel emulation with Rosetta.

.. code-block:: bash

   softwareupdate --install-rosetta


You now need to use your terminal application with Rosetta. Depending on the Mac OS version, your terminal app may be ``iTerm`` or ``Terminal``. ``iTerm`` is found in the Applications folder, while ``Terminal`` may be under the Applications > Utilities directory. Make a copy of this application and name the new copy as "iTerm_Rosetta" or "Terminal_Rosetta". Right click on the renamed terminal app icon, choose "Get Info", and select the "Open using Rosetta" box. To test that you are indeed using the Intel emulator, type the following in your iTerm_Rosetta/Terminal_Rosetta:

.. code-block:: bash

   arch

This command should return ``i386``.

Next, install homebrew with Rosetta:

.. code-block:: bash

   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

With the Intel emulation, Homebrew should be installed at ``/usr/local/bin/brew``. Add the following to your ``.bash_profile``

.. code-block:: bash

   alias brew_i386="/usr/local/bin/brew"

In the future, you will use `brew_i386` as an alternative of `brew` with the Intel emulation.

For completeness only, you might also install Homebrew in your M1/M2/M3 terminal, which should be then installed at ``/opt/homebrew/bin/brew``. Add the following to your ``.bash_profile``

.. code-block:: bash

   alias brew="/opt/homebrew/bin/brew"

Now we will install ``miniconda3`` in Rosetta, but before that, we will have to modify ``.bash_profile`` so we could handle the ``conda`` between M1/M2/M3 and Rosetta separately. Here I assume you already installed ``anaconda`` or ``miniconda`` in your M1/M2/M3 terminal, so the following block should be in your ``.bash_profile``:

.. code-block:: bash

   # >>> conda initialize >>>
   # !! Contents within this block are managed by 'conda init' !!
   __conda_setup="$('/Users/xxxx/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
   if [ $? -eq 0 ]; then
      eval "$__conda_setup"
   else
      if [ -f "/Users/xxxx/anaconda3/etc/profile.d/conda.sh" ]; then
          . "/Users/xxxx/anaconda3/etc/profile.d/conda.sh"
      else
          export PATH="/Users/xxxx/anaconda3/bin:$PATH"
      fi
  fi
  unset __conda_setup
  # <<< conda initialize <<<

Note that the "xxxx" here should be your username. Copy these few lines and paste them into a separate file ``.init_conda_arm64.sh`` in the home directory. Then remove these lines from the file ``.bash_profile``. We will come back to handle the ``.init_conda_arm64.sh`` file later.

Now let's install ``miniconda3`` in Rosetta. First, type the following line in iTerm_Rosetta:

.. code-block:: bash

   curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh > Miniconda3-latest-MacOSX-x86_64.sh

Then type the following and follow instructions to proceed with the installation:

.. code-block:: bash

   bash Miniconda3-latest-MacOSX-x86_64.sh

Once the installation succeed, you will see that the following several new lines have been added to ``.bash_profile``:

.. code-block:: bash

   # >>> conda initialize >>>
   # !! Contents within this block are managed by 'conda init' !!
   __conda_setup="$('/Users/xxxx/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
   if [ $? -eq 0 ]; then
       eval "$__conda_setup"
   else
       if [ -f "/Users/xxxx/miniconda3/etc/profile.d/conda.sh" ]; then
           . "/Users/xxxx/miniconda3/etc/profile.d/conda.sh"
       else
           export PATH="/Users/xxxx/miniconda3/bin:$PATH"
       fi
   fi
   unset __conda_setup
   # <<< conda initialize <<<

Let's cut (i.e., copy and remove) these few lines again and paste them into a separate file ``.init_conda_x86_64.sh`` in the home directory. In the same iTerm_Rosetta, type the following:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

Okay, now we are ready to go ahead and modify ``.bash_profile`` to handle two versions of ``conda`` between M1/M2/M3 and Rosetta terminals. Add the following lines to your ``.bash_profile``:

.. code-block:: bash

   # <<<<<< Added by TR 20220405 <<
   arch_name="$(uname -m)"

   if [ "${arch_name}" = "x86_64" ]; then
       echo "Running on Rosetta using miniconda3"
       source ~/.init_conda_x86_64.sh
   elif [ "${arch_name}" = "arm64" ]; then
       echo "Running on ARM64 using anaconda"
       source ~/.init_conda_arm64.sh
   else
       echo "Unknown architecture: ${arch_name}"
   fi
   # <<<<<<<< end <<<<<<<

Now, when you open iTerm / iTerm_Rosetta, the terminal will display a message to let you know which  ``conda`` version is being used. If you open native iTerm, the terminal will show this:

.. code-block:: bash

   Last login: Sun Nov 24 12:57:39 on ttys001
   Running on ARM64 using anaconda

And when you open iTerm_Rosetta, the terminal will display the following:

.. code-block:: bash

   Last login: Mon Nov 25 10:04:04 on ttys000
   Running on Rosetta using miniconda3



Then, we install ``gfortran`` in iTerm_Rosetta:

.. code-block:: bash

   brew_i386 install gfortran

Everything is ready now, so we should simply install GASTLI with pip. As explained in the :doc:`installation` section, this can be done by creating a conda/venv environment first, activating it and running ``pip install gastli==0.9.1``.

2) Installation with disk image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`François-Xavier Coudert’s github repository <https://github.com/fxcoudert/gfortran-for-macOS>`_ provides gfortran disk images (.dmg) that allow you to install gfortran via an installation wizard on the Apple Silicon (M1, M2, M3) chips. Follow the installation instructions provided in FX Coudert's repository.

Ensure you download and install the version of gfortran compatible with your operating system (M1, M2 or M3) and architecture (ARM64). After installation, verify the default gfortran binary by running the command ``which gfortran`` in your terminal. This command should return the path to the ARM64 binary you just installed.

If the returned path does not point to the installed ARM64 binary, you may need to update your system's PATH to prioritize the correct binary location.






